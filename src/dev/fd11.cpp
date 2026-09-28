#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../quad.hpp"
#include "../ran.hpp"
#include "../sample_mup.hpp"
#include "../solve_quadratic.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <print>
#include <string>

// FD11 (2011), weak-field vacuum limit of (21)--(24).
// Units: km, keV, gauss. Derivation and validation: summary/fd11.md.
// Per photon I=1, V=+2 Im(Ax conj(Ay)), as in FD11 (39).
namespace fd11 {
constexpr double alpha_em = 1 / 137.035999084;
constexpr double B_QED = 4.414005218e13;
constexpr double hbar_c = 1.973269804e-13; // keV km
constexpr double birefringence = alpha_em / (30 * pi * B_QED * B_QED * hbar_c);
using Stokes = double3; // x=Q, y=U, z=V
struct Options {
  double couple = 1e-3, freeze = 1e-3, tolerance = 1e-5;
  bool enabled = true;
  int charge_sign = -1; // electrons: lower sign in FD11 (33)
};
struct Coeff {
  double kappa, c, s; // k0 Delta n, cos(2 phi_B), sin(2 phi_B)
  double twist = 0;
  double3 omega() const { return {kappa * c, kappa * s, twist}; }
};
// Unit quaternion rotating the Stokes sphere. This is a temporary step
// operator, never a photon's electric-field amplitude.
struct Rotation {
  double w = 1;
  double3 v{0, 0, 0};
  static Rotation exponential(double3 a) {
    double angle = a.length();
    double scale =
        angle < 1e-10 ? 0.5 - angle * angle / 48 : std::sin(angle / 2) / angle;
    return {std::cos(angle / 2), scale * a};
  }
  Stokes apply(Stokes s) const {
    auto t = 2 * cross(v, s);
    return s + w * t + cross(v, t);
  }
};
inline Rotation compose(Rotation after, Rotation before) {
  return {after.w * before.w - dot(after.v, before.v),
          after.w * before.v + before.w * after.v + cross(after.v, before.v)};
}
struct Propagator {
  Rotation rotation;
  double phase = 0; // integral of tr(H)/2; discarded after the step
  Stokes apply(Stokes s) const { return rotation.apply(s); }
  double amplitude_change(Stokes before) const {
    // ||(exp(i phase) U - I) A||, expressed entirely in Stokes and U.
    // Stable at small phase: 1-w = |v|^2/(1+w).
    double w = rotation.w;
    double one_minus_w = w > 0 ? dot(rotation.v, rotation.v) / (1 + w) : 1 - w;
    double sh = std::sin(phase / 2);
    double d2 = 2 * (one_minus_w + 2 * w * sh * sh +
                     std::sin(phase) * dot(rotation.v, before));
    return std::sqrt(std::max(0.0, d2));
  }
};
inline Propagator compose(Propagator after, Propagator before) {
  return {compose(after.rotation, before.rotation), after.phase + before.phase};
}
// Fourth-order commutator-free Magnus method; exact for constant H. Unlike
// explicit RK on Q,U,V it has no stability restriction from a large phase.
template <class CoeffFunc>
Propagator magnus(CoeffFunc const &coeff, double x, double h) {
  constexpr double g = 0.28867513459481288225;
  constexpr double a = 0.53867513459481288225, b = -0.03867513459481288225;
  auto c1 = coeff(x + (0.5 - g) * h), c2 = coeff(x + (0.5 + g) * h);
  auto first = Rotation::exponential(h * (a * c1.omega() + b * c2.omega()));
  auto second = Rotation::exponential(h * (b * c1.omega() + a * c2.omega()));
  return {compose(second, first), (11.0 / 12) * h * (c1.kappa + c2.kappa)};
}
inline double overlap(Stokes s, Coeff c, double mu_prime, int charge_sign) {
  double qb = c.c * s.x + c.s * s.y;
  return std::clamp(
      0.25 * ((1 + mu_prime * mu_prime) + (mu_prime * mu_prime - 1) * qb) +
          0.5 * charge_sign * mu_prime * s.z,
      0.0, 1.0);
}
} // namespace fd11

constexpr double M_star = 1.4;                                    // M_sun
constexpr double R_star = 10;                                     // km
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km; // km
constexpr double B_pole = 1e14;                                   // G

using YVector = std::array<double, 3>;

enum class Polarization { O, E };

struct PhotonEvolution {
  static void geodesic(double, YVector const &y, YVector &dydx) {
    auto [r, psi, alpha] = y;
    double L = std::sqrt(1 - rs / r);
    dydx[0] = L * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * L) * (1 - 3 * rs / (2 * r));
  }
  static double event_absorb(double, YVector const &y) { return y[0] - R_star; }
  static double event_escape(double, YVector const &y) { return y[0] - 10000; }
  static auto compute_knots(Boltzmann const &fb, int n_knots) {
    std::vector<double> knots(n_knots);
    for (int i = 0; i < n_knots; ++i) {
      knots[i] = fb.b_min + (fb.b_max - fb.b_min) * i / (n_knots - 1);
    }
    return knots;
  }

  BField const &bfield;
  Boltzmann const &fb;
  StepperDopr5<3, decltype(geodesic)> stepper;
  Ran ran;

  std::vector<double> knots;
  const int n_scan;
  const double orbit_atol, orbit_rtol, quad_atol, quad_rtol;

  Quad quad;

  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;
  fd11::Options pol_options;
  enum class PolStage { Mode, Integrating, Frozen };
  PolStage pol_stage = PolStage::Mode;
  fd11::Stokes stokes{0, 0, 0};
  double pol_start = std::numeric_limits<double>::infinity();
  double frozen_at = std::numeric_limits<double>::infinity();
  double pol_h = 0;
  struct PolSegment {
    double x, end;
    fd11::Stokes before;
  };
  std::vector<PolSegment> pol_segments;
  struct Stats {
    long starts = 0, freezes = 0, accepted = 0, rejected = 0;
    long coupled_scatterings = 0, escaped_mode = 0, escaped_active = 0;
  } stats;

  PhotonEvolution(BField const &bfield, Boltzmann const &fb, int seed, //
                  int n_knots = 4, int n_scan = 4,                     //
                  double orbit_atol = 1e-10, double orbit_rtol = 1e-10,
                  double quad_atol = 1e-8, double quad_rtol = 1e-8)
      : bfield(bfield), fb(fb), stepper(geodesic, orbit_atol, orbit_rtol),
        ran(seed), knots(compute_knots(fb, n_knots)), n_scan(n_scan),
        orbit_atol(orbit_atol), orbit_rtol(orbit_rtol), quad_atol(quad_atol),
        quad_rtol(quad_rtol) {
    stepper.add_event(&event_absorb);
    stepper.add_event(&event_escape);
  }

  void init(double3 _n, double3 _e1, double3 _e2, double _omega_inf,
            Polarization _pol, double r0, double psi0, double alpha0) {
    n = _n;
    e1 = _e1;
    e2 = _e2;
    omega_inf = _omega_inf;
    pol = _pol;
    // Emission and every scattering create an O/E eigenmode (FD11 15--16).
    pol_stage = PolStage::Mode;
    stokes = {0, 0, 0};
    pol_start = frozen_at = std::numeric_limits<double>::infinity();
    pol_h = 0;
    pol_segments.clear();
    stepper.init(0.0, 1e-3 * R_star, {r0, psi0, alpha0});
  }

  struct Screen {
    double3 r, k, x, y, B;
  };
  Screen screen(YVector const &state) const {
    auto [r, psi, alpha] = state;
    auto rh = std::cos(psi) * e1 + std::sin(psi) * e2;
    auto ph = cross(n, rh);
    auto kh = std::cos(alpha) * rh + std::sin(alpha) * ph;
    // The ray-plane normal and n x k form a parallel-transported screen
    // in Schwarzschild spacetime; x cross y = k. No artificial basis spin.
    auto xh = cross(n, kh);
    double mu = std::clamp(rh.z, -1.0, 1.0), rho = std::hypot(rh.x, rh.y);
    auto th = rho > 0 ? double3{rh.x * mu / rho, rh.y * mu / rho, -rho}
                      : double3{std::copysign(1.0, mu), 0, 0};
    auto az = rho > 0 ? double3{-rh.y / rho, rh.x / rho, 0} : double3{0, 1, 0};
    auto bf = bfield.calc_B(r, mu);
    return {rh, kh, xh, n, bf.x * rh + bf.y * th + bf.z * az};
  }
  fd11::Coeff pol_coeff(YVector const &state) const {
    auto f = screen(state);
    double bx = dot(f.B, f.x), by = dot(f.B, f.y), bt2 = bx * bx + by * by;
    double kappa =
        fd11::birefringence * omega_inf / std::sqrt(1 - rs / state[0]) * bt2;
    if (bt2 == 0)
      return {0, 1, 0}; // degenerate eigenmode: deterministic screen
    return {kappa, (bx * bx - by * by) / bt2, 2 * bx * by / bt2};
  }
  fd11::Coeff pol_coeff(double path) const {
    return pol_coeff(stepper.dense_out(path));
  }
  fd11::Stokes eigenmode(fd11::Coeff c) const {
    double sign = pol == Polarization::O ? 1 : -1;
    return {sign * c.c, sign * c.s, 0};
  }
  fd11::Propagator pol_step(double x, double h) const {
    if (h == 0)
      return {};
    auto begin = pol_coeff(x), end = pol_coeff(x + h);
    double angle = std::atan2(begin.c * end.s - begin.s * end.c,
                              begin.c * end.c + begin.s * end.s);
    // In the field-aligned screen the rapid rotation axis is nearly constant:
    // Omega_B=(kappa,0,-(2 phi_B)'). This resolves gradual mode coupling
    // without numerically chasing the much faster dynamical phase in the fixed
    // screen. Large axis turns (near a transverse-field null) use the regular
    // screen.
    if (std::abs(angle) > 0.2)
      return fd11::magnus([&](double s) { return pol_coeff(s); }, x, h);
    auto generator = [&](double path) {
      auto c = pol_coeff(path);
      if (c.kappa == 0)
        return fd11::Coeff{0, 1, 0};
      double dx = 1e-4 * stepper.h_old;
      double a = std::max(stepper.x_old, path - dx);
      double b = std::min(stepper.x_old + stepper.h_old, path + dx);
      auto ca = pol_coeff(a), cb = pol_coeff(b);
      double theta_prime =
          std::atan2(ca.c * cb.s - ca.s * cb.c, ca.c * cb.c + ca.s * cb.s) /
          (b - a);
      return fd11::Coeff{c.kappa, 1, 0, -theta_prime};
    };
    auto p = fd11::magnus(generator, x, h);
    p.rotation =
        fd11::compose(fd11::Rotation::exponential({0, 0, angle}), p.rotation);
    // Conjugate by the starting screen rotation. Avoid half-angle sign flips
    // at atan2's branch cut, which would corrupt the FD11 amplitude criterion.
    auto v = p.rotation.v;
    p.rotation.v = {begin.c * v.x - begin.s * v.y,
                    begin.s * v.x + begin.c * v.y, v.z};
    return p;
  }
  fd11::Stokes polarization_at(double path) const {
    if (path < pol_start)
      return eigenmode(pol_coeff(path));
    if (path >= frozen_at)
      return stokes;
    auto it =
        std::upper_bound(pol_segments.begin(), pol_segments.end(), path,
                         [](double s, PolSegment const &p) { return s < p.x; });
    if (it == pol_segments.begin())
      return eigenmode(pol_coeff(path));
    --it;
    return pol_step(it->x, std::max(0.0, path - it->x)).apply(it->before);
  }
  void advance_polarization(double xl, double xr) {
    pol_segments.clear();
    if (!pol_options.enabled || pol_stage == PolStage::Frozen)
      return;
    double x = xl;
    if (pol_stage == PolStage::Mode) {
      // Locate (34), including photons emitted/scattered beyond the surface.
      auto threshold = [&](double s) {
        return pol_coeff(s).kappa * stepper.dense_out(s)[0] *
                   pol_options.couple -
               1;
      };
      double left = xl, fl = threshold(left);
      bool found = fl <= 0;
      for (int j = 1; !found && j <= n_scan; ++j) {
        double right = xl + (xr - xl) * j / n_scan;
        double fr = threshold(right);
        if (fr <= 0) {
          x = zriddr(threshold, left, right, orbit_atol);
          found = true;
        }
        left = right;
      }
      if (!found)
        return;
      pol_start = x;
      stokes = eigenmode(pol_coeff(x));
      pol_stage = PolStage::Integrating;
      ++stats.starts;
      pol_h = xr - x;
    }
    while (x < xr) {
      double h = std::min(pol_h > 0 ? pol_h : xr - x, xr - x);
      if (x + h == x)
        throw std::runtime_error("Polarization step underflow");
      auto whole = pol_step(x, h);
      auto first = pol_step(x, h / 2), second = pol_step(x + h / 2, h / 2);
      auto half = fd11::compose(second, first);
      auto trial = half.apply(stokes);
      double error = (trial + (-1.0) * whole.apply(stokes)).length() / 15;
      double factor =
          error == 0
              ? 3
              : std::clamp(0.9 * std::pow(pol_options.tolerance / error, 0.2),
                           0.2, 3.0);
      pol_h = h * factor;
      if (error > pol_options.tolerance) {
        ++stats.rejected;
        continue;
      }
      ++stats.accepted;
      pol_segments.push_back({x, x + h / 2, stokes});
      pol_segments.push_back({x + h / 2, x + h, first.apply(stokes)});
      auto end = stepper.dense_out(x + h);
      double change = half.amplitude_change(stokes) * end[0] / h;
      stokes = trial;
      x += h;
      // FD11 (35), including its common phase. The weak-phase guard avoids
      // false freezing when a large phase happens to wrap through 2 pi.
      if (change < pol_options.freeze &&
          pol_coeff(x).kappa * end[0] < 1 && std::cos(end[2]) > 0) {
        pol_stage = PolStage::Frozen;
        frozen_at = x;
        ++stats.freezes;
        return;
      }
    }
  }

  void init_random_radial(double omega_inf, Polarization pol) {
    double3 e1 = ran.point_on_unit_sphere();
    double3 n = ran.unit_perp_to(e1);
    double3 e2 = cross(n, e1);
    init(n, e1, e2, omega_inf, pol, R_star, 0., 0.);
  }

  struct Geometry {
    double x, mu, D, pref, basepref;
  };

  Geometry geo(YVector const &y) const {
    auto [r, psi, alpha] = y;
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = std::clamp(r_hat.z, -1.0, 1.0);
    double3 B_vec = bfield.calc_B(r, muz);
    double B = B_vec.length();
    double3 b = B_vec / B;
    double x = B_to_omega * B * std::sqrt(1 - rs / r) / omega_inf;
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * r_hat.z / rho,
                                          r_hat.y * r_hat.z / rho, -rho}
                                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu = std::clamp(
        b.x * std::cos(alpha) +
            std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat)),
        -1.0, 1.0);
    double D = std::fma(x, x, (mu - 1) * (mu + 1));
    double overlap = (pol == Polarization::E) ? (0.5) : (0.5 * D / (x * x));
    double twist = bfield.Bphi_over_Btheta(muz);
    double pref = (bfield.p + 1) * pi * twist * fb.inv_abs_b_mean * x * x *
                  overlap / (r * std::sqrt(D));
    if (!std::isfinite(pref) || pref < 0) {
      pref = 0;
    }
    double basepref = (bfield.p + 1) * pi * twist * fb.inv_abs_b_mean * x * x /
                      (r * std::sqrt(D));
    if (!std::isfinite(basepref) || basepref < 0)
      basepref = 0;
    return Geometry{
        .x = x,
        .mu = mu,
        .D = D,
        .pref = pref,
        .basepref = basepref,
    };
  }

  // Preserve the original pure-mode opacity exactly. Only the coupling
  // region needs the charge-dependent, beta-dependent overlap of FD11 (33).
  auto polarized_rates(std::array<double, 2> const &betas, Geometry const &g,
                       double path) const {
    auto ret = rate(betas);
    if (!pol_options.enabled || path < pol_start)
      return ret;
    auto s = polarization_at(path);
    auto c = pol_coeff(path);
    for (int i = 0; i < 2; ++i) {
      double mup =
          std::clamp((g.mu - betas[i]) / (1 - betas[i] * g.mu), -1.0, 1.0);
      ret[i] *= fd11::overlap(s, c, mup, pol_options.charge_sign);
    }
    return ret;
  }

  std::array<double, 2> rate(std::array<double, 2> const &betas) const {
    std::array<double, 2> ret{};
    for (int i = 0; i < 2; ++i) {
      double beta = betas[i];
      double f = fb.f(beta);
      if (f == 0)
        continue;
      double t1 = (1 + beta) * (1 - beta);
      ret[i] = f * t1 * std::sqrt(t1);
    }
    return ret;
  }

  // Opacity = A + B Q_B + C V, including both resonance roots.
  double3 opacity_terms(double path) const {
    auto g = geo(stepper.dense_out(path));
    std::array<double, 2> betas;
    if (g.D <= 0 || !solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu,
                                     (1 + g.x) * (1 - g.x), betas))
      return {0, 0, 0};
    auto rates = rate(betas);
    double3 terms{0, 0, 0};
    for (int i = 0; i < 2; ++i) {
      double mu =
          std::clamp((g.mu - betas[i]) / (1 - betas[i] * g.mu), -1.0, 1.0);
      terms = terms + (g.basepref * rates[i]) *
                          double3{0.25 * (1 + mu * mu), 0.25 * (mu * mu - 1),
                                  0.5 * pol_options.charge_sign * mu};
    }
    return terms;
  }
  // In B's rotating screen U_B'=-kappa V-theta' Q_B, theta=2 phi_B.
  // Integrate C V by parts. This is an identity, not phase averaging:
  // int C V = -[C U_B/kappa] + int [(C/kappa)' U_B - C theta' Q_B/kappa].
  // The remaining oscillatory amplitude is suppressed by 1/kappa.
  double fast_optical_depth(double left, double right) {
    double inset = std::min(0.01 * (right - left),
                            1e-8 * (right - left) +
                                16 * std::numeric_limits<double>::epsilon() *
                                    std::max(1.0, right));
    auto boundary = [&](double x, double interior) {
      auto c = pol_coeff(x);
      auto s = polarization_at(x);
      double ub = -c.s * s.x + c.c * s.y;
      // The Boltzmann distribution has a finite jump at beta=0. Quadrature
      // cuts already isolate it; endpoint coefficients must use the limit
      // from inside this interval, irrespective of the root's rounding sign.
      return opacity_terms(interior).z * ub / c.kappa;
    };
    auto integrand = [&](double x) {
      auto c = pol_coeff(x);
      auto s = polarization_at(x);
      double qb = c.c * s.x + c.s * s.y, ub = -c.s * s.x + c.c * s.y;
      auto terms = opacity_terms(x);
      // Centered finite differences of slowly varying coefficients, with a
      // one-sided stencil only at the ends of the orbit's dense domain.
      double dx = std::max(1e-7, 1e-4 * (right - left));
      double a = std::max(left + inset, x - dx);
      double b = std::min(right - inset, x + dx);
      auto ca = pol_coeff(a), cb = pol_coeff(b);
      double ratio_prime =
          (opacity_terms(b).z / cb.kappa - opacity_terms(a).z / ca.kappa) /
          (b - a);
      double theta_prime =
          std::atan2(ca.c * cb.s - ca.s * cb.c, ca.c * cb.c + ca.s * cb.s) /
          (b - a);
      return terms.x + terms.y * qb + ratio_prime * ub -
             terms.z * theta_prime * qb / c.kappa;
    };
    double result = quad.qags(integrand, left, right, quad_atol, quad_rtol) +
                    boundary(left, left + inset) -
                    boundary(right, right - inset);
    if (result < -quad_atol)
      throw std::runtime_error(
          std::format("Negative oscillatory optical depth: value={} "
                      "interval=[{},{}] kappa={} terms={}",
                      result, left, right, pol_coeff(left).kappa,
                      opacity_terms(left).x));
    return std::max(0.0, result);
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
    double Q, U, V;
  } escaped_data;

  struct ScatteredData {
    YVector r_psi_alpha;
    double beta;
  } scattered_data;

  EvolveResult evolve_geodesic() {
    double tau = std::log(ran.U());
    std::vector<double> cuts;

    while (true) {
      stepper.do_step(0.1 * stepper.y_old[0]);
      int event_id = stepper.detect_event();

      double xl = stepper.x_old, xr = stepper.x_old + stepper.h_old;
      advance_polarization(xl, xr);
      cuts.clear(), cuts.push_back(xl), cuts.push_back(xr);
      if (pol_start > xl && pol_start < xr)
        cuts.push_back(pol_start);
      if (frozen_at > xl && frozen_at < xr)
        cuts.push_back(frozen_at);
      double l = xl;
      auto gl = geo(stepper.dense_out(l));
      bool resonant = gl.D > 0;
      for (int i = 1; i <= n_scan; ++i) {
        double r = xl + (xr - xl) * i / n_scan;
        auto gr = geo(stepper.dense_out(r));
        resonant = resonant || gr.D > 0;
        if (gl.D * gr.D <= 0) {
          cuts.push_back(zriddr(
              [&](double x) { return geo(stepper.dense_out(x)).D; }, l, r, 0));
        }
        for (double beta : knots) {
          double ginv = std::sqrt((1 - beta) * (1 + beta));
          auto bound = [beta, ginv](double x, double mu) {
            return x * ginv + beta * mu - 1;
          };
          if (bound(gl.x, gl.mu) * bound(gr.x, gr.mu) <= 0) {
            cuts.push_back(zriddr(
                [&](double x) {
                  auto g = geo(stepper.dense_out(x));
                  return bound(g.x, g.mu);
                },
                l, r, 0));
          }
        }
        l = r;
        gl = gr;
      }
      // Most coupling steps are outside the resonance surface. Only insert
      // dense-polarization joins into quadrature when opacity can be nonzero.
      if (resonant)
        for (auto const &segment : pol_segments)
          if (segment.end > xl && segment.end < xr)
            cuts.push_back(segment.end);
      std::sort(cuts.begin(), cuts.end());
      cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
      for (size_t i = 0; i + 1 < cuts.size(); ++i) {
        double l = cuts[i], r = cuts[i + 1];
        if (geo(stepper.dense_out(std::midpoint(l, r))).D <= 0)
          continue;
        auto integrand = [&](double path_length) {
          auto g = geo(stepper.dense_out(path_length));
          std::array<double, 2> betas;
          if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu,
                               (1 + g.x) * (1 - g.x), betas))
            return 0.0;
          auto rates = polarized_rates(betas, g, path_length);
          double pref = pol_options.enabled && path_length >= pol_start
                            ? g.basepref
                            : g.pref;
          return (rates[0] + rates[1]) * pref;
        };
        bool fast =
            pol_options.enabled && l >= pol_start && r <= frozen_at &&
            pol_coeff(std::midpoint(l, r)).kappa * (r - l) > 256 &&
            geo(stepper.dense_out(l)).D > 1e-8 &&
            geo(stepper.dense_out(r)).D > 1e-8;
        auto optical_depth = [&](double a, double b) {
          return fast && pol_coeff(std::midpoint(a, b)).kappa * (b - a) > 256
                     ? fast_optical_depth(a, b)
                     : quad.qags(integrand, a, b, quad_atol, quad_rtol);
        };
        double dtau = optical_depth(l, r);

        if (tau < 0 && tau + dtau >= 0) {
          double scattered_point = r;
          if (tau + dtau > 0) {
            auto func = [&](double r) {
              if (l == r)
                return tau;
              return tau + optical_depth(l, r);
            };
            auto const &dfunc = integrand;
            scattered_point = rtsafe(func, dfunc, l, r, orbit_atol);
          }

          auto r_psi_alpha = stepper.dense_out(scattered_point);
          auto g = geo(r_psi_alpha);
          std::array<double, 2> betas;
          if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu,
                               (1 + g.x) * (1 - g.x), betas)) {
            throw std::runtime_error("No valid beta at scattering point");
          }
          auto rates = polarized_rates(betas, g, scattered_point);
          double sum_rates = rates[0] + rates[1];
          if (sum_rates <= 0.0 || !std::isfinite(sum_rates)) {
            throw std::runtime_error("Invalid scattering rates");
          }

          scattered_data.r_psi_alpha = r_psi_alpha;
          scattered_data.beta =
              (ran.U() * sum_rates < rates[0]) ? betas[0] : betas[1];
          if (scattered_point >= pol_start)
            ++stats.coupled_scatterings;
          return EvolveResult::Scattered;
        }
        tau += dtau;
      }

      if (event_id == 0) {
        return EvolveResult::Absorbed;
      } else if (event_id == 1) {
        auto [r, psi, alpha] = stepper.y_new;
        double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
        double3 psi_hat = cross(n, r_hat);
        double3 k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
        escaped_data.omega_inf = omega_inf;
        escaped_data.muk = std::clamp(k_hat.z, -1.0, 1.0);
        if (pol_stage == PolStage::Mode) {
          stokes = eigenmode(pol_coeff(stepper.y_new));
          ++stats.escaped_mode;
        } else if (pol_stage == PolStage::Integrating) {
          ++stats.escaped_active;
        }
        // Output screen: x = projected magnetic axis, y = k cross x.
        // This makes photons with different ray planes directly comparable.
        auto f = screen(stepper.y_new);
        double3 observer_x = double3{0, 0, 1} + (-k_hat.z) * k_hat;
        double len = observer_x.length();
        observer_x = len > 1e-12 ? observer_x / len : f.x;
        double c = dot(observer_x, f.x), s = dot(observer_x, f.y);
        double c2 = c * c - s * s, s2 = 2 * c * s;
        escaped_data.Q = c2 * stokes.x + s2 * stokes.y;
        escaped_data.U = -s2 * stokes.x + c2 * stokes.y;
        escaped_data.V = stokes.z;
        return EvolveResult::Escaped;
      }
      stepper.update_old();
    }
  }

  void perform_scattering() {
    auto [r, psi, alpha] = scattered_data.r_psi_alpha;
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double3 B_vec = bfield.calc_B(r, std::clamp(r_hat.z, -1.0, 1.0));
    double B = B_vec.length();
    double3 b = B_vec / B;
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * r_hat.z / rho,
                                          r_hat.y * r_hat.z / rho, -rho}
                                : double3{std::copysign(1.0, r_hat.z), 0, 0};
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu = std::clamp(
        b.x * std::cos(alpha) +
            std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat)),
        -1.0, 1.0);

    double beta = scattered_data.beta;
    double mup_out = sample_mup(ran);
    double mu_out =
        std::clamp((mup_out + beta) / (1 + beta * mup_out), -1.0, 1.0);
    double3 b_cartesian = b.x * r_hat + b.y * theta_hat + b.z * phi_hat;
    double3 t_hat = ran.unit_perp_to(b_cartesian);
    double3 k_out =
        mu_out * b_cartesian + std::sqrt(1 - mu_out * mu_out) * t_hat;
    double cross_len = cross(r_hat, k_out).length();
    double alpha_out = std::atan2(cross_len, dot(k_out, r_hat));
    double3 n_out = (cross_len < 1e-10) ? ran.unit_perp_to(r_hat)
                                        : cross(r_hat, k_out) / cross_len;
    double3 e1_out = r_hat;
    double3 e2_out = cross(n_out, e1_out);
    double omega_inf_out = omega_inf * (1 - beta * mu) / (1 - beta * mu_out);
    Polarization pol_out = (ran.U() < 1 / (1 + mup_out * mup_out))
                               ? Polarization::E
                               : Polarization::O;

    if (!std::isfinite(omega_inf_out)) {
      throw std::runtime_error(
          "Non-finite omega_inf encountered after scattering");
    }

    init(n_out, e1_out, e2_out, omega_inf_out, pol_out, r, 0, alpha_out);
  }
};

BField bfield("table/bfield_t10.txt", B_pole, R_star);
Boltzmann fb(-0.75);

#ifndef FD11_NO_MAIN
int main(int argc, char **argv) try {
  int N = 10000, seed = 7774;
  double energy = 1;
  Polarization initial_mode = Polarization::E;
  fd11::Options options;
  std::string output = "output/fd11.txt";
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto value = [&]() -> std::string {
      if (++i >= argc)
        throw std::invalid_argument("Missing value for " + arg);
      return argv[i];
    };
    if (arg == "--photons")
      N = std::stoi(value());
    else if (arg == "--seed")
      seed = std::stoi(value());
    else if (arg == "--energy")
      energy = std::stod(value());
    else if (arg == "--output")
      output = value();
    else if (arg == "--pol-tol")
      options.tolerance = std::stod(value());
    else if (arg == "--couple")
      options.couple = std::stod(value());
    else if (arg == "--freeze")
      options.freeze = std::stod(value());
    else if (arg == "--no-polarization")
      options.enabled = false;
    else if (arg == "--charge") {
      auto charge = value();
      if (charge != "electron" && charge != "positron")
        throw std::invalid_argument("--charge: electron or positron");
      options.charge_sign = charge == "electron" ? -1 : 1;
    } else if (arg == "--mode") {
      auto mode = value();
      if (mode != "O" && mode != "E")
        throw std::invalid_argument("--mode: O or E");
      initial_mode = mode == "O" ? Polarization::O : Polarization::E;
    } else if (arg == "--help") {
      std::println(
          "FD11 polarization transfer (run from repository root)\n"
          "  --photons N (10000) --seed N (7774) --energy keV (1)\n"
          "  --mode E|O (E) --output output/fd11.txt\n"
          "  --pol-tol 1e-5 --couple 1e-3 --freeze 1e-3 (0 disables freezing)\n"
          "  --charge electron|positron (electron) --no-polarization\n"
          "Output: omega_inf_keV mu_k Q U V; I=1.\n"
          "Screen x=projected magnetic axis, y=k cross x; V=+2 Im(Ax Ay*).");
      return 0;
    } else
      throw std::invalid_argument("Unknown option: " + arg);
  }
  if (N <= 0 || !(energy > 0 && std::isfinite(energy)) ||
      !(options.couple > 0 && std::isfinite(options.couple)) ||
      !(options.tolerance > 0 && std::isfinite(options.tolerance)) ||
      !(options.freeze >= 0 && std::isfinite(options.freeze)))
    throw std::invalid_argument(
        "Invalid count, energy, or polarization tolerance");
  std::filesystem::create_directories("output");
  PhotonEvolution pe(bfield, fb, seed);
  pe.pol_options = options;
  FILE *f = std::fopen(output.c_str(), "w");
  if (!f)
    throw std::runtime_error("Cannot open " + output);
  std::println(f, "# omega_inf_keV mu_k Q U V (I=1; x=projected magnetic axis; "
                  "V=+2Im(AxAy*))");
  long escaped = 0, absorbed = 0, scattered = 0;
  auto start_time = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    pe.init_random_radial(energy, initial_mode);
    while (true) {
      auto result = pe.evolve_geodesic();
      if (result == PhotonEvolution::EvolveResult::Escaped) {
        auto [omega_inf, muk, Q, U, V] = pe.escaped_data;
        std::println(f, "{} {} {} {} {}", omega_inf, muk, Q, U, V);
        ++escaped;
        break;
      } else if (result == PhotonEvolution::EvolveResult::Absorbed) {
        ++absorbed;
        break;
      } else if (result == PhotonEvolution::EvolveResult::Scattered) {
        pe.perform_scattering();
        ++scattered;
      }
    }
  }
  std::fclose(f);
  double seconds = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - start_time)
                       .count();
  std::println(
      stderr, "photons={} escaped={} absorbed={} scatterings={} seconds={:.6f}",
      N, escaped, absorbed, scattered, seconds);
  auto const &s = pe.stats;
  std::println(stderr,
               "pol_starts={} freezes={} accepted={} rejected={} "
               "coupled_scatterings={} escaped_mode={} escaped_active={}",
               s.starts, s.freezes, s.accepted, s.rejected,
               s.coupled_scatterings, s.escaped_mode, s.escaped_active);
  return 0;
} catch (std::exception const &e) {
  std::println(stderr, "fd11: {}", e.what());
  return 1;
}
#endif
