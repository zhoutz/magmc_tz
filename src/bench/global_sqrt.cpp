// Ported from a907a14c98e03c8183a514791c2303a302d15bbe.
// Only Boltzmann API and automatic dense-output preparation are adapted.
#include "benchmark_io.hpp"
#include <format>
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
  Boltzmann b(b0, 4);
  std::vector<double> knots;
  for (int j = 0; j <= 10; ++j)
    knots.push_back(b0 < 0 ? -1 + j * 0.1 : j * 0.1);
  return {[b](double beta) { return b.f(beta); }, b.b_mean, knots};
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
} // namespace od
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

namespace global_sqrt_bench {
double total_optical_depth(optical_bench::Input const &c, std::string const &mode,
                           optical_bench::Counts &work) {
  if (mode != "default" && mode != "matched" && mode != "audit")
    throw std::invalid_argument("Unknown global_sqrt configuration");
  od::Ray ray(c.muz,c.alpha,c.az,c.energy,c.polarization(),c.radius);
  auto dist=od::thermal(c.b0);
  od::Settings cfg;
  cfg.cap=.5; cfg.orbit_tol=1e-10; cfg.tol=1e-8; cfg.scan=1;
  if(mode=="audit") { cfg.cap*=.25; cfg.orbit_tol*=.25; cfg.tol*=.25; }
  od::Counts counts;
  auto record=[&] {work={counts.geometry,counts.density,counts.steps,counts.quadrature,counts.warnings};};
  try { double tau=od::global_integrate(ray,dist,cfg,counts); record(); return tau; }
  catch(...) { record(); throw; }
}
}
#ifndef OPTICAL_BENCH_NO_MAIN
int main(int argc,char **argv) {
  return optical_bench::run(argc,argv,"global_sqrt",global_sqrt_bench::total_optical_depth);
}
#endif
