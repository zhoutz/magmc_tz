#pragma once
// Benchmark-only implementations. Existing production headers are unmodified.
#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
#include "../solve_quadratic.hpp"
#include <chrono>
#include <cstdio>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <memory>
#include <string>

namespace od {
constexpr double radius = 10, rs = 1.4 * schwarzschild_radius_of_sun_in_km;
inline BField field("table/bfield_t10.txt", 1e14, radius);
using State = std::array<double, 3>;
struct Counts {
  size_t geometry = 0, density = 0, steps = 0, quadrature = 0, warnings = 0;
};
// Contract: knots include support boundaries, discontinuities and enough points
// to resolve narrow features. No finite black-box sampler can discover an
// arbitrarily narrow, undeclared peak. f is a density in beta, not momentum.
struct Distribution {
  std::function<double(double)> f;
  double mean;
  std::vector<double> knots;
};
inline Distribution thermal(double b0) {
  Boltzmann b(b0);
  std::vector<double> knots;
  for (int j = 0; j <= 10; ++j)
    knots.push_back(b0 < 0 ? -1 + j * 0.1 : j * 0.1);
  return {[b](double beta) { return b.f(beta); }, b.b_bar(), knots};
}
inline Distribution box(double lo, double hi) {
  return {[=](double b) { return b > lo && b < hi ? 1 / (hi - lo) : 0; },
          (lo + hi) / 2,
          {lo, (lo + hi) / 2, hi}};
}
inline Distribution mixture() {
  auto a = box(-.71, -.69), b = box(-.201, -.199);
  return {[=](double v) { return .65 * a.f(v) + .35 * b.f(v); },
          .65 * a.mean + .35 * b.mean,
          {-.71, -.70, -.69, -.201, -.200, -.199}};
}
struct Ray {
  double muz, alpha, az, energy, r0 = radius, rend = 10000;
  Polarization pol;
  double3 e1, e2, n;
  Ray(double m, double a, double z, double e, Polarization p, double r = radius, double end = 10000)
      : muz(m), alpha(a), az(z), energy(e), r0(r), rend(end), pol(p) {
    e1 = {std::sqrt(1 - m * m), 0, m};
    double3 theta{m, 0, -e1.x}, phi{0, 1, 0};
    e2 = std::cos(z) * theta + std::sin(z) * phi;
    n = cross(e1, e2);
  }
};
inline void orbit(double, State const &y, State &dy) {
  double r = y[0], a = y[2], L = std::sqrt(1 - rs / r);
  dy = {L * std::cos(a), std::sin(a) / r, -std::sin(a) / (r * L) * (1 - 1.5 * rs / r)};
}
struct Geometry {
  double r, muz, x, mu, D, pref;
};
inline Geometry geometry(Ray const &ray, State const &y, Counts &count) {
  ++count.geometry;
  auto [r, psi, alpha] = y;
  auto rh = std::cos(psi) * ray.e1 + std::sin(psi) * ray.e2;
  double muz = std::clamp(rh.z, -1., 1.);
  auto bv = field.calc_B(r, muz);
  double B = bv.length();
  auto b = bv / B;
  double rho = std::hypot(rh.x, rh.y);
  double3 th = rho > 0 ? double3{rh.x * muz / rho, rh.y * muz / rho, -rho}
                       : double3{std::copysign(1., muz), 0, 0};
  double3 ph = rho > 0 ? double3{-rh.y / rho, rh.x / rho, 0} : double3{0, 1, 0};
  double mu =
      b.x * std::cos(alpha) + std::sin(alpha) * (b.y * dot(ray.n, ph) - b.z * dot(ray.n, th));
  mu = std::clamp(mu, -1., 1.);
  double x = B_to_omega * B * std::sqrt(1 - rs / r) / ray.energy;
  return {r,
          muz,
          x,
          mu,
          std::fma(x, x, (mu - 1) * (mu + 1)),
          (field.p + 1) * pi * field.Bphi_over_Btheta(muz) / r};
}
inline double density(Geometry const &g, Distribution const &f, Polarization pol, Counts &c,
                      bool stable) {
  ++c.density;
  if (g.D <= 0) return 0;
  std::array<double, 2> roots;
  if (stable) {
    // Specialized stable quadratic; D is evaluated before coefficient scaling.
    double a = g.x * g.x + g.mu * g.mu;
    double q = g.mu + std::copysign(g.x * std::sqrt(g.D), g.mu);
    roots = {q / a, (1 - g.x) * (1 + g.x) / q};
  } else if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu, (1 - g.x) * (1 + g.x), roots))
    return 0;
  double ret = 0;
  for (double beta : roots) {
    double fval = f.f(beta);
    if (fval == 0) continue;
    double w = 1 - beta * g.mu, muprime = (g.mu - beta) / w;
    double P = pol == Polarization::E ? .5 : .5 * muprime * muprime;
    ret += fval * P * w * w * (1 - beta) * (1 + beta) / std::abs(g.mu - beta);
  }
  ret *= g.pref / std::abs(f.mean);
  if (!std::isfinite(ret)) throw std::runtime_error("nonfinite opacity");
  return ret;
}
struct Settings {
  double cap = .5, tol = 1e-6, orbit_tol = 1e-9;
  double end_length = std::numeric_limits<double>::infinity();
  int scan = 4;
  bool transformed = false, reference = false;
};
// Reuse one workspace per ray; errors are checked, never silently discarded.
struct Quad {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2048);
  ~Quad() { gsl_integration_workspace_free(w); }
  template <class F> double integrate(F &fn, double lo, double hi, double tol, Counts &c) {
    gsl_function gf;
    gf.function = +[](double x, void *p) { return (*static_cast<F *>(p))(x); };
    gf.params = &fn;
    double value, error;
    int status = gsl_integration_qags(&gf, lo, hi, tol, tol, 2048, w, &value, &error);
    ++c.quadrature;
    if (status != GSL_SUCCESS) {
      ++c.warnings;
      // QAGS' extrapolation sometimes encounters the piecewise-linear field
      // table's derivative jumps. Retry without extrapolation before accepting
      // any roundoff-limited estimate. Every retry remains visible in the log.
      if (status != GSL_EROUND || error > 10 * tol * (1 + std::abs(value))) {
        double alternative, alternative_error;
        int other = gsl_integration_qag(&gf, lo, hi, tol, tol, 2048, GSL_INTEG_GAUSS61, w,
                                        &alternative, &alternative_error);
        ++c.quadrature;
        if (other == GSL_SUCCESS || alternative_error < error) {
          status = other;
          value = alternative;
          error = alternative_error;
        }
      }
      // Roundoff at root/table interpolation boundaries can limit QAGS. Accept
      // only an explicit small error estimate, log every such occurrence.
      if (status != GSL_SUCCESS &&
          (status != GSL_EROUND || error > 10 * tol * (1 + std::abs(value))))
        throw std::runtime_error(std::format("quadrature: {} value={:.17g} error={:.17g}",
                                             gsl_strerror(status), value, error));
    }
    if (!std::isfinite(value)) throw std::runtime_error("nonfinite quadrature");
    return value;
  }
};
// Event function H(beta)=x sqrt(1-beta^2)+beta mu-1 has no spurious
// roots from squaring. D=0 marks branch coalescence, independently of f.
inline double boundary(Geometry const &g, double beta) {
  return g.x * std::sqrt((1 - beta) * (1 + beta)) + beta * g.mu - 1;
}
// Differentiate the prescribed piecewise-linear field and the dense orbit.
// This avoids subtracting two nearly equal values of D at a narrow caustic.
struct Dual {
  double v, d;
  Dual(double v = 0, double d = 0) : v(v), d(d) {}
};
inline Dual operator+(Dual a, Dual b) { return {a.v + b.v, a.d + b.d}; }
inline Dual operator-(Dual a, Dual b) { return {a.v - b.v, a.d - b.d}; }
inline Dual operator-(Dual a) { return {-a.v, -a.d}; }
inline Dual operator*(Dual a, Dual b) { return {a.v * b.v, a.d * b.v + a.v * b.d}; }
inline Dual operator/(Dual a, Dual b) { return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)}; }
inline Dual sin(Dual a) { return {std::sin(a.v), std::cos(a.v) * a.d}; }
inline Dual cos(Dual a) { return {std::cos(a.v), -std::sin(a.v) * a.d}; }
inline Dual sqrt(Dual a) {
  double v = std::sqrt(a.v);
  return {v, a.d / (2 * v)};
}
inline Dual pow(Dual a, double p) { return {std::pow(a.v, p), p * std::pow(a.v, p - 1) * a.d}; }
struct Dense {
  double left, right, h;
  State c1, c2, c3, c4, c5;
  template <class T> std::array<T, 3> at(T l) const {
    T s = (l - left) / h, t = 1 - s;
    std::array<T, 3> y;
    for (int i = 0; i < 3; ++i)
      y[i] = c1[i] + s * (c2[i] + t * (c3[i] + s * (c4[i] + t * c5[i])));
    return y;
  }
};
inline double discriminant_derivative(Ray const &ray, Dense const &panel, double l) {
  auto y = panel.at(Dual{l, 1});
  auto r = y[0], psi = y[1], alpha = y[2];
  auto rx = cos(psi) * ray.e1.x + sin(psi) * ray.e2.x;
  auto ry = cos(psi) * ray.e1.y + sin(psi) * ray.e2.y;
  auto rz = cos(psi) * ray.e1.z + sin(psi) * ray.e2.z;
  auto m = rz.v < 0 ? -rz : rz;
  double t = (m.v - field.mu_min) / (field.mu_max - field.mu_min) * (field.mu_num - 1);
  int i = std::clamp(int(t), 0, field.mu_num - 2);
  auto frac = (m - field.mu_min) * ((field.mu_num - 1) / (field.mu_max - field.mu_min)) - i;
  auto fv = (1 - frac) * field.f[i] + frac * field.f[i + 1];
  auto fp = (1 - frac) * field.fp[i] + frac * field.fp[i + 1];
  if (rz.v < 0) fp = -fp;
  auto sth = sqrt(1 - rz * rz), rho = sqrt(rx * rx + ry * ry);
  auto br = -fp, bt = field.p * fv / sth, bp = field.A * pow(fv, 1 / field.p) * bt;
  auto norm = sqrt(br * br + bt * bt + bp * bp);
  auto ndph = (-ray.n.x * ry + ray.n.y * rx) / rho;
  auto ndth = (ray.n.x * rx * rz + ray.n.y * ry * rz) / rho - ray.n.z * rho;
  auto mu = (br * cos(alpha) + sin(alpha) * (bt * ndph - bp * ndth)) / norm;
  auto x = (.5 * field.B_pole * B_to_omega / ray.energy) * pow(radius / r, 2 + field.p) * norm *
           sqrt(1 - rs / r);
  return 2 * (x.v * x.d + mu.v * mu.d);
}
template <class GetGeo, class GetPanel>
inline double smooth_caustic(Ray const &ray, Distribution const &f, Counts &c, GetGeo &geo,
                             GetPanel &panel, double edge, double other, Quad &quad, double tol,
                             double near = 0) {
  double dir = other > edge ? 1. : -1., width = std::abs(other - edge);
  auto ge = geo(edge);
  auto smooth = [&](double u) {
    double delta = u * u, l = edge + dir * delta;
    auto g = geo(l);
    double K;
    if (delta < 1e-5 * ge.r) {
      // Two-point Gauss integration of D' yields a stable divided difference.
      constexpr double v = .28867513459481287;
      double l1 = edge + dir * delta * (.5 - v), l2 = edge + dir * delta * (.5 + v);
      K = dir * .5 *
          (discriminant_derivative(ray, panel(l1), l1) +
           discriminant_derivative(ray, panel(l2), l2));
    } else
      K = (g.D - ge.D) / delta;
    if (K <= 0) throw std::runtime_error("non-simple caustic: split at stationary D first");
    double D = delta * K, sq = g.x * std::sqrt(D), den = g.x * g.x + g.mu * g.mu;
    // Use consistent D for the second root near coalescence.
    std::array<double, 2> roots{(g.mu - sq) / den, (g.mu + sq) / den};
    double ret = 0;
    ++c.density;
    for (double beta : roots) {
      double w = 1 - beta * g.mu, P = ray.pol == Polarization::E ? .5 : .5 * D / (g.x * g.x);
      ret += f.f(beta) * P * w * (1 - beta) * (1 + beta) * g.x;
    }
    // |mu-beta| = sqrt(D)*(1-beta*mu)/x; 2u/sqrt(D)=2/sqrt(K).
    return ret * g.pref / std::abs(f.mean) * 2 / std::sqrt(K);
  };
  return quad.integrate(smooth, std::sqrt(near), std::sqrt(width), tol, c);
}
struct Cut {
  double s;
  bool singular;
};
template <class F>
inline void find_events(std::vector<Cut> &cuts, double a, double b, double fa, double fm, double fb,
                        F fn, bool singular) {
  double m = (a + b) / 2;
  auto bracket = [&](double x, double y, double fx, double fy) {
    if (fx * fy < 0) cuts.push_back({zriddr(fn, x, y, 1e-14), singular});
    if (fx == 0) cuts.push_back({x, singular});
    if (fy == 0) cuts.push_back({y, singular});
  };
  bracket(a, m, fa, fm);
  bracket(m, b, fm, fb);
  // An extremum can enclose two crossings even when all sampled values have
  // the same sign. Search it explicitly, instead of trusting endpoint signs.
  // Unimodality within this small scan cell is a resolution assumption; the
  // cap/scan convergence checks and grazing-ray regression test exercise it.
  double curvature = fa - 2 * fm + fb;
  bool vertex_inside = curvature != 0 && std::abs(fa - fb) < 2 * std::abs(curvature);
  if (fa * fm > 0 && fm * fb > 0 && ((fm - fa) * (fb - fm) < 0 || vertex_inside)) {
    double sign = curvature > 0 ? 1. : -1.;
    constexpr double golden = .6180339887498948482;
    double lo = a, hi = b, x = hi - golden * (hi - lo), y = lo + golden * (hi - lo);
    double fx = sign * fn(x), fy = sign * fn(y);
    for (int it = 0; it < 48; ++it) {
      if (fx < fy) {
        hi = y;
        y = x;
        fy = fx;
        x = hi - golden * (hi - lo);
        fx = sign * fn(x);
      } else {
        lo = x;
        x = y;
        fx = fy;
        y = lo + golden * (hi - lo);
        fy = sign * fn(y);
      }
    }
    double e = (lo + hi) / 2, fe = fn(e);
    bracket(a, e, fa, fe);
    bracket(e, b, fe, fb);
  }
}
inline double integrate(Ray const &ray, Distribution const &f, Settings cfg, Counts &c) {
  StepperDopr5<3, decltype(orbit)> step(orbit, cfg.orbit_tol, cfg.orbit_tol);
  step.init(0, .01, {ray.r0, 0, ray.alpha});
  Quad quad;
  double tau = 0;
  for (size_t it = 0; it < 100000; ++it) {
    step.do_step(cfg.cap * step.y_old[0]);
    step.prepare_dense();
    ++c.steps;
    double left = step.x_old, right = left + step.h_old;
    Dense dense{left,        right,       step.h_dense, step.rcont1,
                step.rcont2, step.rcont3, step.rcont4,  step.rcont5};
    auto state = [&](double s) { return step.dense_out(s); };
    auto geo = [&](double s) { return geometry(ray, state(s), c); };
    bool done = false;
    if (right >= cfg.end_length) {
      right = cfg.end_length;
      done = true;
    }
    if (step.y_new[0] >= ray.rend || step.y_new[0] <= radius) {
      double end = step.y_new[0] >= ray.rend ? ray.rend : radius;
      double hit =
          zriddr([&](double s) { return state(s)[0] - end; }, left, left + step.h_old, 1e-14);
      right = std::min(right, hit);
      done = true;
    }
    std::vector<Cut> cuts{{left, false}, {right, false}};
    // Several samples per orbit step detect exit and re-entry; the cap and
    // scan density must be checked for convergence for a new field geometry.
    double a = left;
    auto ga = geo(a);
    for (int j = 1; j <= cfg.scan; ++j) {
      double b = left + (right - left) * j / cfg.scan;
      auto gb = geo(b), gm = geo((a + b) / 2);
      find_events(cuts, a, b, ga.D, gm.D, gb.D, [&](double s) { return geo(s).D; }, true);
      for (double beta : f.knots)
        if (std::abs(beta) < 1) {
          find_events(
              cuts, a, b, boundary(ga, beta), boundary(gm, beta), boundary(gb, beta),
              [&](double s) { return boundary(geo(s), beta); }, false);
        }
      a = b;
      ga = gb;
    }
    std::sort(cuts.begin(), cuts.end(), [](auto a, auto b) { return a.s < b.s; });
    std::vector<Cut> unique;
    for (auto v : cuts) {
      if (!unique.empty() &&
          v.s - unique.back().s < 4 * std::numeric_limits<double>::epsilon() * (1 + std::abs(v.s)))
        unique.back().singular |= v.singular;
      else
        unique.push_back(v);
    }
    for (size_t j = 1; j < unique.size(); ++j) {
      auto a = unique[j - 1], b = unique[j];
      if (b.s <= a.s) continue;
      auto mid = geo((a.s + b.s) / 2);
      if (cfg.transformed && mid.D <= 0) continue;
      // Reference deliberately retains zero panels, as ref.cpp does.
      if (!cfg.reference && density(mid, f, ray.pol, c, true) == 0) {
        // This skip is valid only for intervals bounded by every support edge;
        // interior declared features are split as well. Do not infer support
        // from an underflowing f: check kinematic roots against knot span.
        if (mid.D <= 0) continue;
        double sq = mid.x * std::sqrt(mid.D), den = mid.x * mid.x + mid.mu * mid.mu;
        double v1 = (mid.mu - sq) / den, v2 = (mid.mu + sq) / den;
        if ((v1 < f.knots.front() || v1 > f.knots.back()) &&
            (v2 < f.knots.front() || v2 > f.knots.back()))
          continue;
      }
      auto raw = [&](double s) { return density(geo(s), f, ray.pol, c, !cfg.reference); };
      double nearest = std::numeric_limits<double>::infinity(), distance = nearest;
      for (auto cut : unique)
        if (cut.singular) {
          double d = std::min(std::abs(cut.s - a.s), std::abs(cut.s - b.s));
          if (d < distance) {
            distance = d;
            nearest = cut.s;
          }
        }
      if (cfg.transformed && std::isfinite(nearest)) {
        // l=l_c +/- u^2 removes the E-mode square-root singularity.
        // Both endpoints singular: split at midpoint before transforming.
        auto part = [&](double edge, double other, double near = 0) {
          auto panel = [&](double) -> Dense const & { return dense; };
          return smooth_caustic(ray, f, c, geo, panel, edge, other, quad, cfg.tol, near);
        };
        if (a.singular && b.singular)
          tau += part(a.s, (a.s + b.s) / 2) + part(b.s, (a.s + b.s) / 2);
        else
          tau += part(nearest, nearest <= a.s ? b.s : a.s, distance);
      } else if (!a.singular && !b.singular && distance < b.s - a.s && distance > 0) {
        // A support cut may lie almost on D=0. Raw QAGS extrapolation can
        // mistake this *finite* endpoint for the singularity and double-count
        // its small remaining tail. Geometric panels resolve that separation.
        double dir = nearest <= a.s ? 1. : -1.,
               far = std::max(std::abs(a.s - nearest), std::abs(b.s - nearest));
        for (double d = distance; d < far;) {
          double next = std::min(far, 4 * d);
          double x = nearest + dir * d, y = nearest + dir * next;
          tau += quad.integrate(raw, std::min(x, y), std::max(x, y), cfg.tol, c);
          d = next;
        }
      } else
        tau += quad.integrate(raw, a.s, b.s, cfg.tol, c);
    }
    if (done) return tau;
    step.update_old();
  }
  throw std::runtime_error("orbit step limit");
}
inline double limited(Ray const &ray, Distribution const &f, double cap, double tol, Counts &c) {
  auto deriv = [&](double s, std::array<double, 4> const &y, std::array<double, 4> &dy) {
    State p{y[0], y[1], y[2]}, dp;
    orbit(s, p, dp);
    dy = {dp[0], dp[1], dp[2], density(geometry(ray, p, c), f, ray.pol, c, false)};
  };
  StepperDopr5<4, decltype(deriv)> step(deriv, tol, tol);
  step.add_event([&](double, auto const &y) { return y[0] - ray.rend; });
  step.add_event([&](double, auto const &y) { return y[0] - radius; });
  step.init(0, .01, {ray.r0, 0, ray.alpha, 0});
  for (size_t i = 0; i < 100000; ++i) {
    step.do_step(cap * step.y_old[0]);
    ++c.steps;
    if (step.detect_event() >= 0) return step.y_new[3];
    step.update_old();
  }
  throw std::runtime_error("limited step limit");
}
template <class Integrator> int run(std::string name, Integrator solver, int argc, char **argv) {
  gsl_set_error_handler_off();
  std::string suite = argc > 1 ? argv[1] : "radial";
  std::string suffix = argc > 2 ? argv[2] : "";
  std::string path = "output/" + name + (suite == "radial" ? "" : "_" + suite) + suffix + ".txt";
  FILE *fp = std::fopen(path.c_str(), "w");
  if (!fp) throw std::runtime_error("output open");
  Counts count;
  int id = 0, fail = 0;
  auto start = std::chrono::steady_clock::now();
  auto calc = [&](Ray const &ray, Distribution const &dist) {
    try {
      return solver(ray, dist, count);
    } catch (std::exception const &e) {
      std::fprintf(stderr, "%s case %d: %s\n", name.c_str(), id, e.what());
      ++fail;
      return std::numeric_limits<double>::quiet_NaN();
    }
  };
  if (suite == "radial") {
    for (int k = 1; k <= 9; ++k)
      for (int m = 0; m < 10; ++m)
        for (double e : {.01, .1, 1., 10., 100.})
          for (auto p : {Polarization::E, Polarization::O}) {
            auto dist = thermal(-k * .1);
            Ray ray(m * .1, 0, 0, e, p);
            double tau = calc(ray, dist);
            std::fprintf(fp, "%.2f %.2f %.2f %d %.16e\n", -k * .1, m * .1, e, p == Polarization::E,
                         tau);
            ++id;
          }
  } else if (suite == "nonradial" || suite == "stress") {
    // Northern/southern hemispheres, both polarizations, genuinely 3D orbital
    // planes, near-tangent launch, inward rays that turn or hit the surface.
    for (double m : {-.8, -.3, .0, .4, .85})
      for (double a : {.3, 1.1, 1.55, 2.2})
        for (double az : {.2, 1.4, 3.2})
          for (double e : {.1, 3., 70.})
            for (auto p : {Polarization::E, Polarization::O})
              for (int d = 0; d < 3; ++d) {
                Distribution dist = suite == "nonradial" ? thermal(d == 0   ? -.1
                                                                   : d == 1 ? -.5
                                                                            : -.9)
                                    : d == 0             ? box(-.401, -.399)
                                    : d == 1             ? mixture()
                                                         : box(-.75, .25);
                Ray ray(m, a, az, e, p, a > 1.6 ? 30 : radius);
                double tau = calc(ray, dist);
                std::fprintf(fp, "%d %.8g %.8g %.8g %.8g %.8g %d %d %.16e\n", id, m, a, az, ray.r0,
                             e, p == Polarization::E, d, tau);
                ++id;
              }
  } else
    throw std::runtime_error("suite must be radial, nonradial or stress");
  double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::fclose(fp);
  std::printf("%s %s cases=%d seconds=%.9f geometry=%zu density=%zu steps=%zu quadrature=%zu "
              "roundoff=%zu failures=%d\n",
              name.c_str(), suite.c_str(), id, sec, count.geometry, count.density, count.steps,
              count.quadrature, count.warnings, fail);
  return fail ? 1 : 0;
}
} // namespace od
