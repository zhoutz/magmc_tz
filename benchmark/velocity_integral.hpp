#pragma once
#include "benchmark_common.hpp"

// Independent reference-style beta quadrature with an explicit support
// interval. Radius is obtained by a fresh bracketed solve at each sample.
inline Result velocity_integral(RadialResonance const &ray,double tolerance,
                                 double rmax=10000,double low=-1,double high=0) {
  Result result;
  if(ray.x(R_star)<=1 || rmax<=R_star) return result;
  double a=std::max(low,ray.beta(R_star)), b=std::min(high,ray.beta(rmax));
  if(a>=b) return result;
  auto f=[&](double beta) {
    double G=ray.g(beta);
    // A support endpoint reconstructed through beta -> g can round a few
    // ulps outside the original radius bracket. Use its endpoint limit.
    double r;
    double xlo=ray.x(R_star),xhi=ray.x(rmax);
    constexpr double fuzz=32*std::numeric_limits<double>::epsilon();
    if(G>=xlo && G<=xlo*(1+fuzz))r=R_star;
    else if(G<=xhi && G>=xhi*(1-fuzz))r=rmax;
    else r=zriddr([&](double r){return std::log(ray.x(r)/G);},R_star,rmax,2e-12);
    double m=(ray.mu-beta)/(1-beta*ray.mu), P=ray.pol==Polarization::E ? .5 : .5*m*m;
    return ray.prefactor*ray.fb.f(beta)*P*(1-beta*ray.mu)/(ray.lapse(r)*ray.Q(r));
  };
  std::vector<double> panels={a,b};
  for(double k:{.125,.25,.5,1.,2.,4.,8.}) {
    double v=k*ray.fb.b0;
    if(a<v && v<b)panels.push_back(v);
  }
  for(double v:{-.9,-.99,-.999,-.9999})if(a<v && v<b)panels.push_back(v);
  std::sort(panels.begin(),panels.end());
  // Keep true support endpoints; discard only interior knots separated
  // from their neighbour/endpoints by floating-point roundoff.
  std::vector<double> clean={a};
  for(double v:panels) if(v-clean.back()>32*std::numeric_limits<double>::epsilon() && b-v>32*std::numeric_limits<double>::epsilon())clean.push_back(v);
  clean.push_back(b);panels=std::move(clean);
  for(size_t i=1;i<panels.size();++i) result.tau+=quadrature::adaptive(f,panels[i-1],panels[i],tolerance*.1/(panels.size()-1),tolerance,result.evaluations);
  result.steps=panels.size()-1;
  return result;
}
