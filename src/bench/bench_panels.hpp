// Adapted from origin/codex-run2; benchmark-local implementation.
#pragma once
#include "bench_physics.hpp"
#include "bench_quadrature.hpp"
#include "../roots.hpp"
#include <vector>

namespace transport {
struct Boundary { double t; bool singular=false; };

// A path segment is parameterized by t in [0,1], NOT by radius. Events
// remain valid through radial turning points and changing magnetic angles.
template<class Path, class Material=Medium> struct ResonancePanels {
  Material const &medium;
  Path const &path;
  double length, low, high;
  std::vector<Boundary> boundaries;
  long evaluations=0;
  int caustics=0, double_branches=0;

  ResonancePoint point(double t) const { return medium.at(path(t)); }
  double surface(double t,double beta) const {
    auto p=point(t);
    return beta==2 ? p.logx-.5*std::log(std::max(1e-30,1-p.mu*p.mu)) : p.surface(beta);
  }

  ResonancePanels(Material const &m,Path const &p,double h,std::vector<double> const &velocities,
                  int probes,double lo=-1,double hi=1)
      :medium(m),path(p),length(h),low(lo),high(hi) {
    boundaries={{0,false},{1,false}};
    std::vector<ResonancePoint> samples;
    for(int i=0;i<=probes;++i)samples.push_back(point(double(i)/probes));
    std::vector<double> surfaces=velocities;surfaces.push_back(2); // D=0
    for(double beta:surfaces) {
      std::vector<std::pair<double,double>> nodes;
      for(int i=0;i<=probes;++i) {
        auto p=samples[i];
        nodes.push_back({double(i)/probes,beta==2 ? p.logx-.5*std::log(std::max(1e-30,1-p.mu*p.mu)) : p.surface(beta)});
      }
      // Endpoints can have equal signs around a thin resonance pocket.
      // Locate interior extrema before root bracketing, in addition to
      // sampling. Probe refinement is independently checked in validation.
      for(int i=1;i<probes;++i) {
        double s1=nodes[i].second-nodes[i-1].second,s2=nodes[i+1].second-nodes[i].second;
        if(s1*s2>=0)continue;
        double a=nodes[i-1].first,b=nodes[i+1].first,sgn=s1<0 ? 1. : -1.;
        constexpr double g=.6180339887498948482;
        double c=b-g*(b-a),d=a+g*(b-a),fc=sgn*surface(c,beta),fd=sgn*surface(d,beta);
        for(int k=0;k<48;++k) {
          if(fc<fd){b=d;d=c;fd=fc;c=b-g*(b-a);fc=sgn*surface(c,beta);}
          else {a=c;c=d;fc=fd;d=a+g*(b-a);fd=sgn*surface(d,beta);}
        }
        double t=(a+b)/2;
        nodes.push_back({t,surface(t,beta)});
      }
      std::sort(nodes.begin(),nodes.end());
      for(size_t i=1;i<nodes.size();++i) {
        auto [a,fa]=nodes[i-1];auto [b,fb]=nodes[i];
        if(a==b)continue;
        if(fa==0)boundaries.push_back({a,beta==2});
        if(fb==0)boundaries.push_back({b,beta==2});
        if(std::signbit(fa)==std::signbit(fb))continue;
        double t=zriddr([&](double t){return surface(t,beta);},a,b,3e-15);
        boundaries.push_back({t,beta==2});
      }
    }
    std::sort(boundaries.begin(),boundaries.end(),[](auto a,auto b){return a.t<b.t;});
    std::vector<Boundary> clean;
    for(auto b:boundaries) {
      if(!clean.empty() && b.t-clean.back().t<8e-15)clean.back().singular|=b.singular;
      else clean.push_back(b);
    }
    boundaries=std::move(clean);
    for(auto b:boundaries)caustics+=b.singular;
  }

  double integrate(Boundary left,Boundary right,double atol,double rtol) {
    double a=left.t,b=right.t,w=b-a;
    if(!(w>0))return 0;
    auto mid=point((a+b)/2);double D=mid.discriminant();
    if(D<=0)return 0; // all D=0 crossings have been split
    int contributing=0;
    for(double beta:medium.branches(mid,D))contributing+=beta>std::max(low,medium.distribution.low) && beta<std::min(high,medium.distribution.high);
    if(contributing==0)return 0;
    double_branches+=contributing==2;

    // Local Taylor coefficients anchor D at a located simple caustic.
    // They avoid subtracting nearly equal O(1) numbers in tiny neighborhoods
    // after the square-root change of variable. No angular specialization.
    auto coefficients=[&](double t,double direction) {
      auto y=path(t),p=medium.at(y);auto rates=medium.rates(y);
      double first=2*std::exp(2*p.logx)*rates[0]+2*p.mu*rates[1];
      double dl=1e-4*y[0];auto dy=geodesic(y);State yp=y,ym=y;
      for(int k=0;k<3;++k){yp[k]+=dl*dy[k];ym[k]-=dl*dy[k];}
      auto pp=medium.at(yp),pm=medium.at(ym);
      auto rp=medium.rates(yp),rm=medium.rates(ym);
      double second=(2*std::exp(2*pp.logx)*rp[0]+2*pp.mu*rp[1]
                    -2*std::exp(2*pm.logx)*rm[0]-2*pm.mu*rm[1])/(2*dl);
      return std::array<double,2>{direction*first*length,second*length*length};
    };
    std::array<double,2> lc{},rc{};
    if(left.singular)lc=coefficients(a,1);
    if(right.singular)rc=coefficients(b,-1);
    auto integrand=[&](double z) {
      // sin^2 maps both possible endpoint singularities at once.
      double sa=std::sin(.5*pi*z),sb=std::cos(.5*pi*z);
      double da=w*sa*sa,db=w*sb*sb;
      double t=z<.5 ? a+da : b-db;
      auto p=point(t);double D=p.discriminant();
      // A panel may be extremely narrow in D. Scaling the Taylor region by
      // panel width alone leaves most of such a panel dominated by roundoff.
      if(left.singular && da*length<1e-5*path(a)[0] && std::abs(lc[0]*da)<1e-5)
        D=lc[0]*da+.5*lc[1]*da*da;
      if(right.singular && db*length<1e-5*path(b)[0] && std::abs(rc[0]*db)<1e-5)
        D=rc[0]*db+.5*rc[1]*db*db;
      return medium.opacity(p,D,low,high)*length*w*pi*sa*sb;
    };
    return quadrature::adaptive(integrand,0,1,atol,rtol,evaluations);
  }
};
}
