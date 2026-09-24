// Adapted from origin/codex-run2; benchmark-local implementation.
#pragma once
#include "bench_distribution.hpp"
#include <algorithm>

namespace transport {
using State = std::array<double,3>;

inline State geodesic(State const &y) {
  double r=y[0], alpha=y[2], L=std::sqrt(1-rs/r);
  return {L*std::cos(alpha),std::sin(alpha)/r,
          -std::sin(alpha)/(r*L)*(1-1.5*rs/r)};
}

struct ResonancePoint {
  double logx, mu, prefactor;
  double discriminant() const { return std::expm1(2*logx)+mu*mu; }
  double surface(double beta) const {
    return logx-std::log1p(-beta*mu)+.5*std::log1p(-beta*beta);
  }
};

struct Medium {
  BField const &field;
  Distribution const &distribution;
  Photon const &photon;
  mutable long evaluations=0;
  mutable long opacity_evaluations=0;

  ResonancePoint at(State const &y) const {
    ++evaluations;
    auto [r,psi,alpha]=y;
    double3 rhat=std::cos(psi)*photon.e1+std::sin(psi)*photon.e2;
    // Unit-vector roundoff must not place a pole outside the magnetic table.
    double muz=std::clamp(rhat.z,-1.,1.);
    auto B=field.calc_B(r,muz);double norm=B.length();auto b=B/norm;
    double rho=std::hypot(rhat.x,rhat.y);
    double3 theta=rho>1e-14 ? double3{rhat.x*muz/rho,rhat.y*muz/rho,-rho} : double3{std::copysign(1.,muz),0,0};
    double3 phi=rho>1e-14 ? double3{-rhat.y/rho,rhat.x/rho,0} : double3{0,1,0};
    double mu=b.x*std::cos(alpha)+std::sin(alpha)*(b.y*dot(photon.n,phi)-b.z*dot(photon.n,theta));
    return {std::log(norm*B_to_omega/photon.omega_inf)+.5*std::log1p(-rs/r),
            std::clamp(mu,-1.,1.),
            (field.p+1)*pi*field.Bphi_over_Btheta(muz)/(r*std::abs(distribution.b_bar()))};
  }

  std::array<double,2> branches(ResonancePoint const &p,double D) const {
    if(!(D>0))return {NAN,NAN};
    double x=std::exp(p.logx), den=1+D, s=x*std::sqrt(D);
    // Stable quadratic roots: the smaller numerator is recovered by Vieta.
    double large=(p.mu+std::copysign(s,p.mu))/den;
    // Use the same D as the caustic regularization, including in Vieta's
    // product. Otherwise one root retains the cancelled, unregularized D.
    double small=(p.mu*p.mu-D)/(den*large);
    if(!std::isfinite(small))small=(p.mu-std::copysign(s,p.mu))/den;
    return {std::min(large,small),std::max(large,small)};
  }

  double opacity(ResonancePoint const &p,double D,double low=-1,double high=1) const {
    ++opacity_evaluations;
    if(!(D>0) || p.prefactor==0)return 0;
    double x=std::exp(p.logx),s=std::sqrt(D), sum=0;
    // |mu-beta| = sqrt(D)*(1-beta*mu)/x. This form is stable
    // when both roots approach mu; both contributing branches are retained.
    double overlap=photon.pol==Polarization::E ? .5 : .5*D/(x*x);
    for(double b:branches(p,D)) {
      if(b<low || b>high)continue;
      double f=distribution.f(b);
      if(f==0)continue;
      sum+=f*overlap*(1-b*p.mu)*(1-b*b)*x/s;
    }
    double value=p.prefactor*sum;
    if(!std::isfinite(value) || value<0)throw std::runtime_error("Invalid resonant opacity");
    return value;
  }

  // Directional derivatives along the full geodesic, including variation of
  // B, mu, photon direction and gravitational frequency. No radial formula.
  std::array<double,2> rates(State const &y) const {
    auto dy=geodesic(y);double h=1e-5*y[0];State a=y,b=y;
    for(int i=0;i<3;++i){a[i]-=h*dy[i];b[i]+=h*dy[i];}
    auto pa=at(a),pb=at(b);
    return {(pb.logx-pa.logx)/(2*h),(pb.mu-pa.mu)/(2*h)};
  }
};

inline Photon make_photon(double r,double muz,double alpha,double azimuth,double energy,Polarization pol) {
  double sth=std::sqrt(std::max(0.,1-muz*muz));
  double3 er{sth,0,muz},etheta{muz,0,-sth},ephi{0,1,0};
  double3 tangent=std::cos(azimuth)*etheta+std::sin(azimuth)*ephi;
  return {cross(er,tangent),er,tangent,r,0,alpha,energy,pol};
}
}
