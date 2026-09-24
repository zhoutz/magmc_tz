#include "benchmark_common.hpp"
#include "../src/transport_fast.hpp"

// Convergence audit of the same historical algorithm: no tail clipping,
// all magnetic interpolation knots, a smaller geometry step.
int main(int argc,char **argv) {
  Config defaults;defaults.tol=1e-9;defaults.parameter=.025;
  return benchmark_main(argc,argv,"transport_fast_audit",[](auto const &field,auto const &c,auto const &cfg,auto const &) {
    Boltzmann fb(c.b0);
    auto p=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
    fast_transport::TransportOptions options;
    options.quadrature_rtol=cfg.tol;options.quadrature_atol=.01*cfg.tol;
    options.guard_extremum_bounds=true;
    options.geometry_rtol=options.geometry_atol=1e-11;
    options.max_step_fraction=cfg.parameter;
    options.fast_tail_u=std::numeric_limits<double>::infinity();options.fast_split_field_knots=true;
    auto r=fast_transport::propagate(field,fb,p,options);
    return Result{r.tau,static_cast<long>(r.stats.field_evaluations),
                  static_cast<long>(r.stats.geometry_steps),static_cast<long>(15*r.stats.quadrature_intervals),"ok"};
  },defaults," geometry_rtol=1e-11 geometry_atol=1e-11 quadrature_atol=0.01*tolerance initial_step_policy=parameter*r_or_target_bound tail_u=inf split_field_knots=true guard_extremum_bounds=true");
}
