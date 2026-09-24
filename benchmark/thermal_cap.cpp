#include "../output/benchmark_common.hpp"
int main(int argc,char **argv) {
  return benchmark_main(argc,argv,"thermal_cap",[](auto const &f,auto const &c,auto const &cfg,auto edges){return spatial_ode(f,c,cfg,edges,2);});
}
