#pragma once

#include "radial_resonance.hpp"
#include "quadrature.hpp"
#include <vector>

// Spatial integration for outward radial rays with mu >= 0, beta0 < 0.
// r = r_end*(1-u^2) removes the equatorial E-mode square-root singularity.
// Seeded panels make the thermal layer visible BEFORE adaptive quadrature.
// No velocity tail is discarded. This is not a nonradial transport solver.
struct RadialOpticalDepth {
  RadialResonance ray;
  double r_end=R_star, u_surface=0;
  std::vector<double> panels;
  mutable long evaluations=0;

  RadialOpticalDepth(BField const &field, Boltzmann const &fb, double muz,
                     double energy, Polarization pol)
      : ray(field,fb,muz,energy,pol) {
    if (!(fb.b0<0 && ray.mu>=0))
      throw std::invalid_argument("RadialOpticalDepth requires beta0<0 and mu>=0");
    if (!(energy>0 && std::isfinite(energy))) throw std::invalid_argument("Invalid photon energy");
    if (ray.x(R_star)<=1) return;
    // Find the support edge; do not impose the benchmark's 10000 km cutoff.
    double hi=2*R_star;
    while (ray.x(hi)>1) {
      hi*=2;
      if (!std::isfinite(hi)) throw std::runtime_error("Cannot bracket resonance edge");
    }
    r_end=zriddr([&](double r){return std::log(ray.x(r));},R_star,hi,2e-13);
    u_surface=std::sqrt((r_end-R_star)/r_end);
    panels={0,u_surface};
    // Resolve both a cold core and a relativistic tail; nodes select spatial
    // panels only, and are not a cutoff of the integration domain.
    std::vector<double> velocities={-.9,-.99,-.999,-.9999};
    for(double k:{.125,.25,.5,1.,2.,4.,8.}) if(k*std::abs(fb.b0)<1) velocities.push_back(k*fb.b0);
    for(double b:velocities) {
      double G=ray.g(b);
      if(G>=ray.x(R_star)) continue;
      double r=zriddr([&](double r){return std::log(ray.x(r)/G);},R_star,r_end,2e-13);
      panels.push_back(std::sqrt((r_end-r)/r_end));
    }
    std::sort(panels.begin(),panels.end());
    panels.erase(std::unique(panels.begin(),panels.end()),panels.end());
  }

  double density(double u) const {
    double u2=u*u, r=r_end*(1-u2), k=rs/r_end;
    // Compute x^2-1 relative to the known edge, avoiding cancellation when
    // u approaches zero. Evaluating x(r)-1 directly destroys the endpoint.
    double logx=-(ray.q+.5)*std::log1p(-u2)+.5*std::log1p(-u2/(1-k));
    double c=std::expm1(2*logx), X=std::exp(logx);
    double b=-c/(ray.mu+X*std::sqrt(c+ray.mu*ray.mu));
    double f=ray.fb.f(b);
    if(f==0) return 0;
    double m=(ray.mu-b)/(1-b*ray.mu);
    double P=ray.pol==Polarization::E ? .5 : .5*m*m;
    return ray.prefactor/r/ray.lapse(r)*f*P*(1-b*ray.mu)*(1-b*ray.mu)*(1-b*b)*
           (2*r_end*u/std::abs(ray.mu-b));
  }

  // Cumulative optical depth from the surface to r_max, useful for event roots.
  double integrate_to(double r_max=10000, double atol=1e-9, double rtol=1e-8) const {
    if (!(atol>0 && rtol>0)) throw std::invalid_argument("Invalid quadrature tolerance");
    evaluations=0;
    if(r_max<=R_star || u_surface==0) return 0;
    double umin=r_max>=r_end ? 0 : std::sqrt((r_end-r_max)/r_end);
    double total=0;
    for(size_t i=1;i<panels.size();++i) {
      double a=std::max(umin,panels[i-1]), b=panels[i];
      if(a>=b) continue;
      total+=quadrature::adaptive([&](double u){return density(u);},a,b,
                                 atol/(panels.size()-1),rtol,evaluations);
    }
    return total;
  }
};
