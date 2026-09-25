#include "od_global.hpp"
#include <fstream>
// Deterministic supplementary checks: cumulative depth, fixed optical-depth
// draws, and a very thin allowed island around a tangency of D(l).
int main() {
  gsl_set_error_handler_off();
  od::Settings fast, ref;
  fast.scan = 1;
  fast.orbit_tol = 1e-10;
  fast.tol = 1e-8;
  ref.cap = .005;
  ref.orbit_tol = 1e-12;
  ref.tol = 1e-9;
  ref.scan = 1;
  ref.transformed = true;
  ref.reference = true;
  od::Counts counts;
  FILE *out = std::fopen("output/cumulative.txt", "w");
  std::fprintf(out, "# id length tau_global tau_fine abs_difference\n");
  int id = 0;
  double max_diff = 0, max_location = 0;
  FILE *events = std::fopen("output/scattering_events.txt", "w");
  std::fprintf(events,
               "# id U tau_draw length_global length_fine abs_difference (-1 means escape)\n");
  for (double m : {-.7, .2, .8})
    for (double a : {.7, 2.2}) {
      od::Ray r(m, a, 1.4, 3., Polarization::E, a > 2 ? 30 : 10);
      auto f = od::thermal(-.5);
      for (double length : {10., 30., 60., 90., 150., 300.}) {
        auto fs = fast, rs = ref;
        fs.end_length = rs.end_length = length;
        double x = od::global_integrate(r, f, fs, counts), y = od::integrate(r, f, rs, counts);
        max_diff = std::max(max_diff, std::abs(x - y));
        std::fprintf(out, "%d %.16e %.16e %.16e %.16e\n", id, length, x, y, std::abs(x - y));
      }
      for (double U : {.2, .5, .8}) {
        double target = -std::log(U);
        auto locate = [&](bool reference) {
          auto settings = reference ? ref : fast;
          settings.end_length = 300;
          auto value = [&](double length) {
            settings.end_length = length;
            return reference ? od::integrate(r, f, settings, counts)
                             : od::global_integrate(r, f, settings, counts);
          };
          if (value(300) < target) return -1.;
          double lo = 0, hi = 300;
          for (int j = 0; j < 32; ++j) {
            double mid = (lo + hi) / 2;
            if (value(mid) < target)
              lo = mid;
            else
              hi = mid;
          }
          return (lo + hi) / 2;
        };
        double x = locate(false), y = locate(true);
        if ((x < 0) != (y < 0)) throw std::runtime_error("scatter/escape mismatch");
        max_location = std::max(max_location, std::abs(x - y));
        std::fprintf(events, "%d %.8g %.16e %.16e %.16e %.16e\n", id, U, target, x, y,
                     std::abs(x - y));
      }
      ++id;
    }
  std::fclose(out);
  std::fclose(events);
  // Find the maximum Ecrit(l)=B_to_omega B L / sqrt(1-mu^2), independently
  // on a densely stepped geodesic. Choosing E just below its maximum produces
  // two D=0 crossings separated by a very thin D>0 island.
  // Reflect the geometry so the resonating test population has negative mean
  // velocity, consistent with the benchmark's current-carrying electrons.
  od::Ray r(-.4, 2.2, pi + 1.4, 1., Polarization::E, 30.);
  StepperDopr5<3, decltype(od::orbit)> step(od::orbit, 1e-13, 1e-13);
  step.init(0, .001, {r.r0, 0, r.alpha});
  std::vector<od::Dense> path;
  double best = 0, best_l = 0, best_mu = 0;
  while (step.y_old[0] < 1000) {
    step.do_step(.001 * step.y_old[0]);
    step.prepare_dense();
    od::Dense d{step.x_old,   step.x_old + step.h_old,
                step.h_dense, step.rcont1,
                step.rcont2,  step.rcont3,
                step.rcont4,  step.rcont5};
    path.push_back(d);
    auto g = od::geometry(r, step.y_old, counts);
    double crit = g.x / std::sqrt(1 - g.mu * g.mu);
    if (crit > best) {
      best = crit;
      best_l = step.x_old;
      best_mu = g.mu;
    }
    step.update_old();
  }
  auto geom = [&](double l) {
    auto d = std::lower_bound(path.begin(), path.end(), l,
                              [](auto const &p, double v) { return p.right < v; });
    return od::geometry(r, d->at(l), counts);
  };
  double lo = best_l - 1, hi = best_l + 1;
  for (int i = 0; i < 80; ++i) {
    double x = lo + (hi - lo) / 3, y = hi - (hi - lo) / 3;
    auto gx = geom(x), gy = geom(y);
    if (gx.x / std::sqrt(1 - gx.mu * gx.mu) < gy.x / std::sqrt(1 - gy.mu * gy.mu))
      lo = x;
    else
      hi = y;
  }
  best_l = (lo + hi) / 2;
  auto gb = geom(best_l);
  best = gb.x / std::sqrt(1 - gb.mu * gb.mu);
  best_mu = gb.mu;
  auto f = od::box(std::max(-.999, best_mu - .05), std::min(.999, best_mu + .05));
  FILE *grazing = std::fopen("output/grazing.txt", "w");
  std::fprintf(grazing, "# epsilon energy beta_mean tau_global tau_fine rel_difference\n");
  double max_grazing = 0;
  for (double eps : {1e-3, 1e-5, 1e-7}) {
    r.energy = best * (1 - eps);
    auto fs = fast, rs = ref;
    // This nearly tangent geometry is ill-conditioned; tighten orbit accuracy.
    fs.orbit_tol = rs.orbit_tol = 1e-13;
    fs.tol = rs.tol = 1e-6;
    double x = od::global_integrate(r, f, fs, counts), y = od::integrate(r, f, rs, counts);
    if (!(x > 0 && y > 0)) throw std::runtime_error("missed grazing island");
    double err = std::abs(x - y) / y;
    max_grazing = std::max(max_grazing, err);
    std::fprintf(grazing, "%.8g %.16e %.16e %.16e %.16e %.16e\n", eps, r.energy, best_mu, x, y,
                 err);
    std::fflush(grazing);
  }
  std::fclose(grazing);
  std::printf("cumulative_max_abs=%.12g scatter_length_max_abs_km=%.12g grazing_max_rel=%.12g "
              "roundoff=%zu\n",
              max_diff, max_location, max_grazing, counts.warnings);
  return max_diff < 1e-5 && max_location < 1e-3 && max_grazing < 1e-4 ? 0 : 1;
}
