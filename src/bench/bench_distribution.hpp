#pragma once
#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../photon.hpp"
#include "../roots.hpp"
#include "../qags.hpp"
#include <array>
#include <functional>
#include <string>

inline constexpr double R_star=10, rs=1.4*schwarzschild_radius_of_sun_in_km;
inline constexpr double B_pole=1e14;

namespace transport {
// General velocity PDF interface. The solver never inspects a temperature or
// assumes a functional form. Supply ALL discontinuities, support edges and
// landmarks resolving narrow peaks in knots. A finite sampler cannot discover
// an arbitrarily narrow, unadvertised peak in a black-box function.
// bands are {low, median/peak, high} for FT07's controller, one per peak.
struct Distribution {
  std::string name;
  std::function<double(double)> pdf;
  double low, high, mean;
  std::vector<double> knots;
  std::vector<std::array<double,3>> bands;
  std::vector<double> reference_knots; // optional coarser partition for raw QAGS
  double f(double beta) const {
    return beta>low && beta<high ? pdf(beta) : 0.;
  }
  double b_bar() const { return mean; }
};

// Only this benchmark adapter knows about Boltzmann. Quantiles are cached by
// the driver, outside the timed transport loops. No distribution is truncated.
inline Distribution boltzmann_distribution(double b0) {
  Boltzmann fb(b0), positive(std::abs(b0));
  auto quantile=[&](double probability) {
    return zriddr([&](double b) {
      return qags([&](double v){return positive.f(v);},0.,b,1e-13,1e-12)-probability;
    },0.,1.,1e-14);
  };
  double a=quantile(.001), m=quantile(.5), b=quantile(.999);
  std::vector<double> knots={0};
  // The PDF supplies these landmarks; they are not part of the algorithm.
  for(double prob:{.001,.01,.1,.25,.5,.75,.9,.99,.999})
    knots.push_back(std::copysign(quantile(prob),b0));
  for(double fraction:{.1,.25,.5,.75,.9,.99,.999})
    knots.push_back(std::copysign(b+(1-b)*fraction,b0));
  std::array<double,3> band=b0<0 ? std::array<double,3>{-b,-m,-a} : std::array<double,3>{a,m,b};
  return {"boltzmann",[fb](double v){return fb.f(v);},fb.b_min,fb.b_max,fb.b_bar(),knots,{band},{0}};
}

inline Distribution hat_distribution(bool bimodal=false) {
  // A 2e-4-wide peak, or two disconnected peaks with 20% and 80% weights.
  std::vector<std::array<double,3>> parts=bimodal ?
    std::vector<std::array<double,3>>{{-.701,-.699,.2},{-.1201,-.1199,.8}} :
    std::vector<std::array<double,3>>{{-.2001,-.1999,1.}};
  Distribution d{bimodal ? "two_hats" : "narrow_hat",{},-1,0,0,{}, {},{}};
  for(auto [lo,hi,w]:parts) {
    d.mean+=w*(lo+hi)/2;
    d.knots.insert(d.knots.end(),{lo,(lo+hi)/2,hi});
    d.bands.push_back({lo+.001*(hi-lo),(lo+hi)/2,hi-.001*(hi-lo)});
  }
  d.pdf=[parts](double v){double f=0;for(auto [lo,hi,w]:parts)if(v>lo && v<hi)f+=w/(hi-lo);return f;};
  return d;
}

inline Distribution compact_distribution() {
  // Smooth compact, asymmetric mixture on BOTH signs of beta; each component
  // is 15/(16*w)*(1-z^2)^2 for |z|<1. Integral=1, mean=center.
  const std::vector<std::array<double,3>> parts={{-.35,.0002,.8},{.5,.03,.2}};
  Distribution d{"compact_mixture",{},-1,1,-.18,{}, {},{}};
  for(auto [c,w,weight]:parts) {
    for(double z:{-1.,-.5,0.,.5,1.})d.knots.push_back(c+w*z);
    d.bands.push_back({c-w,c,c+w});
  }
  d.pdf=[parts](double v){double f=0;for(auto [c,w,weight]:parts){double z=(v-c)/w;if(std::abs(z)<1)f+=weight*15/(16*w)*std::pow(1-z*z,2);}return f;};
  return d;
}
}
