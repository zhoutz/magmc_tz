#include "../output/benchmark_common.hpp"
#include "../src/radial_optical_depth.hpp"
int main(int argc,char **argv) {
  return benchmark_main(argc,argv,"spatial_regularized",[](auto const &f,auto const &c,auto const &cfg,auto){
    Boltzmann fb(c.b0);
    RadialOpticalDepth integral(f,fb,c.muz,c.energy,c.pol ? Polarization::E : Polarization::O);
    double tau=integral.integrate_to(10000,cfg.tol*.1,cfg.tol);
    return Result{tau,integral.evaluations,static_cast<long>(integral.panels.size())-1,0};
  });
}
