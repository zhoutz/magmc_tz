// FT07's step control with the existing DOPRI5 optical-depth integrator.
// Physics and GR corrections follow doc/resonant.typ, not FT07's flat-space model.
#include "../output/benchmark_common.hpp"
int main(int argc,char **argv) {
  return benchmark_main(argc,argv,"ft07",[](auto const &f,auto const &c,auto const &cfg,auto edges){return spatial_ode(f,c,cfg,edges,1);});
}
