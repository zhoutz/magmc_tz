#include "benchmark_common.hpp"
#include "../src/transport_fast.hpp"

int main(int argc,char **argv) {
  // Match the historical defaults, including the metadata written by the harness.
  Config defaults; defaults.tol=1e-7;
  return benchmark_main(argc,argv,"transport_fast_guarded",[&](auto const &field,auto const &c,auto const &cfg,auto const &) {
    Boltzmann fb(c.b0);
    auto p=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
    fast_transport::TransportOptions options;
    options.guard_extremum_bounds=true;
      options.quadrature_rtol=cfg.tol;
      options.quadrature_atol=.01*cfg.tol;
      options.max_step_fraction=cfg.parameter;
    auto r=fast_transport::propagate(field,fb,p,options);
    return Result{r.tau,static_cast<long>(r.stats.field_evaluations),
                  static_cast<long>(r.stats.geometry_steps),static_cast<long>(15*r.stats.quadrature_intervals),"ok"};
  },defaults," geometry_rtol=1e-9 geometry_atol=1e-10 quadrature_atol=0.01*tolerance initial_step_policy=parameter*r_or_target_bound tail_u=12 split_field_knots=false guard_extremum_bounds=true");
}
