#include "../output/velocity_integral.hpp"
int main(int argc,char **argv) {
  return benchmark_main(argc,argv,"velocity",[](auto const &f,auto const &c,auto const &cfg,auto){
    Boltzmann fb(c.b0);
    RadialResonance ray(f,fb,c.muz,c.energy,c.pol ? Polarization::E : Polarization::O);
    return velocity_integral(ray,cfg.tol);
  });
}
