#pragma once
#include "od_common.hpp"
namespace od {
// Build a reusable dense geodesic, then integrate only between physical
// velocity/caustic events. Numerical orbit steps are not quadrature boundaries.
// rend is a finite endpoint, not an assumed permanent escape at the first D=0.
inline double global_integrate(Ray const &ray, Distribution const &f, Settings cfg, Counts &c) {
  StepperDopr5<3, decltype(orbit)> step(orbit, cfg.orbit_tol, cfg.orbit_tol);
  step.init(0, .01, {ray.r0, 0, ray.alpha});
  std::vector<Dense> path;
  std::vector<Cut> cuts{{0, false}};
  for (size_t it = 0; it < 100000; ++it) {
    step.do_step(cfg.cap * step.y_old[0]);
    step.prepare_dense();
    ++c.steps;
    Dense d{step.x_old,   step.x_old + step.h_old,
            step.h_dense, step.rcont1,
            step.rcont2,  step.rcont3,
            step.rcont4,  step.rcont5};
    bool done = false;
    if (d.right >= cfg.end_length) {
      d.right = cfg.end_length;
      done = true;
    }
    if (step.y_new[0] >= ray.rend || step.y_new[0] <= radius) {
      double end = step.y_new[0] >= ray.rend ? ray.rend : radius;
      double hit = zriddr([&](double s) { return d.at(s)[0] - end; }, d.left,
                          step.x_old + step.h_old, 1e-14);
      d.right = std::min(d.right, hit);
      done = true;
    }
    path.push_back(d);
    auto geo = [&](double s) { return geometry(ray, d.at(s), c); };
    double a = d.left;
    auto ga = geo(a);
    for (int j = 1; j <= cfg.scan; ++j) {
      double b = d.left + (d.right - d.left) * j / cfg.scan;
      auto gb = geo(b), gm = geo((a + b) / 2);
      find_events(cuts, a, b, ga.D, gm.D, gb.D, [&](double s) { return geo(s).D; }, true);
      for (double beta : f.knots)
        if (std::abs(beta) < 1)
          find_events(
              cuts, a, b, boundary(ga, beta), boundary(gm, beta), boundary(gb, beta),
              [&](double s) { return boundary(geo(s), beta); }, false);
      a = b;
      ga = gb;
    }
    if (done) break;
    step.update_old();
    if (it == 99999) throw std::runtime_error("global orbit step limit");
  }
  cuts.push_back({path.back().right, false});
  std::sort(cuts.begin(), cuts.end(), [](auto a, auto b) { return a.s < b.s; });
  std::vector<Cut> unique;
  for (auto v : cuts) {
    if (!unique.empty() &&
        v.s - unique.back().s < 4 * std::numeric_limits<double>::epsilon() * (1 + std::abs(v.s)))
      unique.back().singular |= v.singular;
    else
      unique.push_back(v);
  }
  auto panel = [&](double s) -> Dense const & {
    auto p = std::lower_bound(path.begin(), path.end(), s,
                              [](Dense const &d, double v) { return d.right < v; });
    return p == path.end() ? path.back() : *p;
  };
  auto geo = [&](double s) { return geometry(ray, panel(s).at(s), c); };
  auto raw = [&](double s) { return density(geo(s), f, ray.pol, c, true); };
  double tau = 0;
  Quad quad;
  for (size_t j = 1; j < unique.size(); ++j) {
    auto a = unique[j - 1], b = unique[j];
    if (b.s <= a.s) continue;
    auto g = geo((a.s + b.s) / 2);
    if (g.D <= 0) continue;
    double sq = g.x * std::sqrt(g.D), den = g.x * g.x + g.mu * g.mu;
    double v1 = (g.mu - sq) / den, v2 = (g.mu + sq) / den;
    if ((v1 < f.knots.front() || v1 > f.knots.back()) &&
        (v2 < f.knots.front() || v2 > f.knots.back()))
      continue;
    double nearest = std::numeric_limits<double>::infinity(), distance = nearest;
    for (auto cut : unique)
      if (cut.singular) {
        double d = std::min(std::abs(cut.s - a.s), std::abs(cut.s - b.s));
        if (d < distance) {
          distance = d;
          nearest = cut.s;
        }
      }
    if (std::isfinite(nearest)) {
      auto part = [&](double edge, double other, double near = 0) {
        return smooth_caustic(ray, f, c, geo, panel, edge, other, quad, cfg.tol, near);
      };
      if (a.singular && b.singular)
        tau += part(a.s, (a.s + b.s) / 2) + part(b.s, (a.s + b.s) / 2);
      else
        tau += part(nearest, nearest <= a.s ? b.s : a.s, distance);
    } else
      tau += quad.integrate(raw, a.s, b.s, cfg.tol, c);
  }
  return tau;
}
} // namespace od
