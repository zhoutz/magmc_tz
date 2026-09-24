#pragma once

// Propagate only the three geometrical variables with DOPRI5. Optical depth is
// integrated separately on dense trajectories, split at resonance surfaces.
#include "quad.hpp"
#include "transport_fast_dopr5.hpp"
#include "transport_fast_physics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Port of f373eb507956aeb22eff5e16577573d26b171efc.
// Fast-only propagation with direct Gauss--Kronrod quadrature. The optional
// extremum guard rejects out-of-step roots caused by integer-bound rounding.
namespace fast_transport {
struct TransportOptions {
  double geometry_rtol = 1e-9;
  double geometry_atol = 1e-10;
  double quadrature_rtol = 1e-7;
  double quadrature_atol = 1e-9;
  double max_step_fraction = 0.1;
  // Fast mode omits only the exp(-u^2) Boltzmann tail beyond this landmark.
  // +infinity disables this optimization for convergence checks.
  double fast_tail_u = 12;
  bool fast_split_field_knots = false;
  // Retained historical option for rejecting rounded out-of-step extrema.
  bool guard_extremum_bounds = false;
  std::size_t max_geometry_steps = 200000;
  std::size_t max_quadrature_intervals = 8192;
};

enum class TransportTermination { escaped, absorbed, scattered };

struct TransportStats {
  std::size_t geometry_steps = 0;
  std::size_t geometry_rhs_evaluations = 0;
  std::size_t geometry_evaluations = 0;
  std::size_t field_evaluations = 0;
  std::size_t rate_evaluations = 0;
  std::size_t quadrature_intervals = 0;
  std::size_t event_roots = 0;
  double estimated_quadrature_error = 0;
};

struct TransportResult {
  YVector state{};
  double distance = 0;
  double tau = 0;
  TransportTermination termination = TransportTermination::escaped;
  TransportStats stats{};
};

namespace resonance_transport_detail {
using GeometryState = std::array<double, 3>;

struct GeometryRHS {
  TransportStats &stats;
  void operator()(double, GeometryState const &y, GeometryState &dy) const {
    ++stats.geometry_rhs_evaluations;
    double const r = y[0], alpha = y[2];
    double const f = std::sqrt(1 - rs / r);
    dy[0] = f * std::cos(alpha);
    dy[1] = std::sin(alpha) / r;
    dy[2] = -std::sin(alpha) / (r * f) * (1 - 1.5 * rs / r);
  }
};

struct Surface {
  enum Kind { beta, discriminant, radius, magnetic_knot } kind;
  double value;
};

struct Cut {
  double s;
};
} // namespace resonance_transport_detail

// tau_target is an absolute accumulated optical depth; initial[3] is retained.
// No RNG is used. A caller may pass -log(U) and scatter at the returned state.
inline TransportResult transport_fast(PhotonEvolution const &photon, YVector initial,
                                      TransportOptions const &options = {},
                                      double tau_target = std::numeric_limits<double>::infinity(),
                                      double r_escape = 1000 * R_star) {
  using namespace resonance_transport_detail;
  if (!(options.fast_tail_u >= 8))
    throw std::invalid_argument("fast_tail_u must be at least 8 (or +infinity)");
  if (!(options.geometry_rtol > 0) || !(options.geometry_atol > 0) ||
      !(options.quadrature_rtol > 0) || !(options.quadrature_atol > 0) ||
      !(options.max_step_fraction > 0 && options.max_step_fraction <= .25) ||
      !std::isfinite(initial[3]) || std::isnan(tau_target) || !(r_escape > R_star) ||
      initial[0] < R_star * (1 - 32 * std::numeric_limits<double>::epsilon()))
    throw std::invalid_argument("invalid resonance transport initial state or options");
  TransportResult result;
  result.state = initial;
  result.tau = initial[3];
  if (result.tau >= tau_target) {
    result.termination = TransportTermination::scattered;
    return result;
  }
  if (initial[0] >= r_escape) return result;
  if (initial[0] <= R_star && std::cos(initial[2]) < 0) {
    result.termination = TransportTermination::absorbed;
    return result;
  }
  GeometryRHS rhs{result.stats};
  StepperDopr5<3, GeometryRHS> stepper(rhs, options.geometry_atol, options.geometry_rtol);
  ++result.stats.field_evaluations;
  auto const radial_geometry = photon.geometry(initial, true);
  // Restart at the local geometry scale. Resonance splitting already protects
  // the optical-depth layer; two tiny startup steps after every scatter add no
  // accuracy. The geometry RK error controller can reject this first proposal.
  double first_step = options.max_step_fraction * initial[0];
  if (std::isfinite(tau_target)) {
    double const initial_rate = photon.rate(radial_geometry);
    if (initial_rate > 0)
      first_step = std::min(first_step, 1.5 * (tau_target - initial[3]) / initial_rate);
  }
  stepper.init(0, first_step, {initial[0], initial[1], initial[2]});

  // Thermal landmarks and an explicit, configurable Boltzmann-tail boundary.
  std::vector<Surface> surfaces{{Surface::discriminant, 0},
                                {Surface::beta, 0},
                                {Surface::radius, R_star},
                                {Surface::radius, r_escape}};
  for (double t : {1., 4., std::isfinite(options.fast_tail_u) ? options.fast_tail_u : 8.}) {
    double const gm1 = t * t / photon.fb.a;
    double const beta = std::copysign(std::sqrt(gm1 * (gm1 + 2)) / (1 + gm1), photon.fb.b0);
    if (std::abs(beta) < 1) surfaces.push_back({Surface::beta, beta});
  }

  // A radial ray keeps its angular field and direction fixed: do not rebuild
  // a spherical basis or interpolate the same field at every quadrature node.
  bool const radial = initial[2] == 0 || initial[2] == pi;
  double const tail_gamma = std::isfinite(options.fast_tail_u)
                                ? 1 + options.fast_tail_u * options.fast_tail_u / photon.fb.a
                                : 1;
  double const tail_beta =
      std::copysign(std::sqrt((tail_gamma - 1) * (tail_gamma + 1)) / tail_gamma, photon.fb.b0);
  auto in_population = [&](ResonanceGeometry const &g) {
    if (!(g.discriminant > 0)) return false;
    if (!std::isfinite(options.fast_tail_u)) return true;
    double const edge = tail_gamma * (1 - tail_beta * g.mu);
    double lower = std::min(1., edge), upper = std::max(1., edge);
    if (tail_beta * g.mu > 0 && std::abs(g.mu) < std::abs(tail_beta))
      lower = std::sqrt((1 - g.mu) * (1 + g.mu));
    return g.x >= lower && g.x <= upper;
  };
  auto locate = [&](auto const &fn, double a, double b) {
    return zriddr(fn, a, b, 2 * std::numeric_limits<double>::epsilon());
  };
  auto rate = [&](ResonanceGeometry const &g) {
    ++result.stats.rate_evaluations;
    std::array<double, 2> betas{};
    auto weights = photon.resonance_weights(g, photon.pol, betas);
    return (weights[0] + weights[1]) * (photon.bfield.p + 1) * pi * g.current /
           (std::abs(photon.fb.b_bar()) * g.r);
  };

  for (std::size_t count = 0; count < options.max_geometry_steps; ++count) {
    stepper.h_old = std::min(stepper.h_old, options.max_step_fraction * stepper.y_old[0]);
    stepper.do_step();
    stepper.prepare_dense();
    ++result.stats.geometry_steps;
    double const h = stepper.h_old;
    auto state = [&](double s) -> YVector {
      // Evaluate the dense polynomial directly: computing x_old+s*h would lose
      // the small offsets needed to resolve an endpoint resonance singularity.
      GeometryState y;
      double const s1 = 1 - s;
      for (int j = 0; j < 3; ++j)
        y[j] = stepper.rcont1[j] +
               s * (stepper.rcont2[j] +
                    s1 * (stepper.rcont3[j] + s * (stepper.rcont4[j] + s1 * stepper.rcont5[j])));
      return {y[0], y[1], y[2], result.tau};
    };
    auto geometry = [&](double s) {
      ++result.stats.geometry_evaluations;
      if (radial) {
        auto g = radial_geometry;
        g.r = state(s)[0];
        g.x *= std::pow(initial[0] / g.r, 2 + photon.bfield.p) *
               std::sqrt((1 - rs / g.r) / (1 - rs / initial[0]));
        g.discriminant = std::fma(g.x, g.x, (g.mu - 1) * (g.mu + 1));
        return g;
      }
      ++result.stats.field_evaluations;
      return photon.geometry(state(s), true);
    };
    constexpr int probes = 2;
    std::array<ResonanceGeometry, probes + 1> sample{}, minus{}, plus{};
    constexpr double derivative_delta = 1e-5;
    for (int i = 0; i <= probes; ++i) {
      double const s = double(i) / probes;
      sample[i] = geometry(s);
      minus[i] = i == 0 ? sample[i] : geometry(std::max(0., s - derivative_delta));
      plus[i] = i == probes ? sample[i] : geometry(std::min(1., s + derivative_delta));
    }
    std::vector<Cut> cuts{{0}, {1}};
    double stop = 1;
    bool stopped = false;
    TransportTermination stop_reason = TransportTermination::escaped;
    auto add_cut = [&](double s, Surface const &surface) {
      ++result.stats.event_roots;
      cuts.push_back({s});
      if (surface.kind == Surface::radius && s > 16 * std::numeric_limits<double>::epsilon() &&
          s <= stop) {
        stop = s;
        stopped = true;
        stop_reason = surface.value == R_star ? TransportTermination::absorbed
                                              : TransportTermination::escaped;
      }
    };
    bool possible =
        std::any_of(sample.begin(), sample.end(), [](auto const &g) { return g.discriminant > 0; });
    auto find_surface = [&](Surface const &surface) {
      double const beta_sqrt =
          surface.kind == Surface::beta ? std::sqrt((1 - surface.value) * (1 + surface.value)) : 0;
      auto value = [&](ResonanceGeometry const &g) {
        if (surface.kind == Surface::discriminant && g.discriminant > 0) possible = true;
        if (surface.kind == Surface::beta)
          return std::fma(g.x, beta_sqrt, surface.value * g.mu - 1);
        switch (surface.kind) {
        case Surface::discriminant:
          return g.discriminant;
        case Surface::radius:
          return g.r - surface.value;
        case Surface::magnetic_knot:
          return g.muz - surface.value;
        case Surface::beta:
          break; // returned above using cached beta_sqrt
        }
        return 0.;
      };
      auto fn = [&](double s) {
        if (surface.kind == Surface::radius) return state(s)[0] - surface.value;
        if (surface.kind == Surface::magnetic_knot) {
          double psi = state(s)[1];
          return photon.e1.z * std::cos(psi) + photon.e2.z * std::sin(psi) - surface.value;
        }
        return value(geometry(s));
      };
      for (int i = 0; i < probes; ++i) {
        double const a = double(i) / probes, b = double(i + 1) / probes;
        double const fa = value(sample[i]), fb = value(sample[i + 1]);
        double const da = value(plus[i]) - value(minus[i]);
        double const db = value(plus[i + 1]) - value(minus[i + 1]);
        std::array<double, 3> local{a, b};
        std::size_t local_size = 2;
        // Endpoint signs alone miss double crossings. Search a stationary point
        // whenever its derivative is bracketed by the geometrical probes.
        if (da * db < 0) {
          bool const minimum = da < 0;
          double lo = a, hi = b;
          constexpr double ratio = .6180339887498948482;
          double c = hi - ratio * (hi - lo), d = lo + ratio * (hi - lo);
          double fc = fn(c), fd = fn(d);
          for (int k = 0; k < 55 && hi - lo > 4e-15; ++k) {
            if ((fc < fd) == minimum) {
              hi = d;
              d = c;
              fd = fc;
              c = hi - ratio * (hi - lo);
              fc = fn(c);
            } else {
              lo = c;
              c = d;
              fc = fd;
              d = lo + ratio * (hi - lo);
              fd = fn(d);
            }
          }
          double const extreme = .5 * (lo + hi);
          local = {a, extreme, b};
          local_size = 3;
          // Even without a zero this point may be a narrow near-tangent peak.
          cuts.push_back({extreme});
        }
        for (std::size_t j = 1; j < local_size; ++j) {
          double const l = local[j - 1], r = local[j];
          double const fl = l == a ? fa : fn(l), fr = r == b ? fb : fn(r);
          if (fl == 0) add_cut(l, surface);
          if (fr == 0) add_cut(r, surface);
          if (std::signbit(fl) != std::signbit(fr)) add_cut(locate(fn, l, r), surface);
        }
      }
    };
    find_surface(surfaces[0]); // D, including interior extrema
    find_surface(surfaces[2]); // stellar surface
    find_surface(surfaces[3]); // escape surface
    if (!possible) {
      result.distance = stepper.x_old + stop * h;
      result.state = state(stop);
      result.state[3] = result.tau;
      if (stopped) {
        result.termination = stop_reason;
        return result;
      }
      stepper.update_old();
      continue;
    }
    for (auto const &surface : surfaces)
      if (surface.kind == Surface::beta) find_surface(surface);

    std::sort(cuts.begin(), cuts.end(), [](Cut a, Cut b) { return a.s < b.s; });
    bool active = false;
    for (std::size_t j = 1; j < cuts.size(); ++j) {
      double a = cuts[j - 1].s, b = std::min(stop, cuts[j].s);
      if (b > a && in_population(geometry(.5 * (a + b)))) {
        active = true;
        break;
      }
    }
    if (!active) {
      result.distance = stepper.x_old + stop * h;
      result.state = state(stop);
      result.state[3] = result.tau;
      if (stopped) {
        result.termination = stop_reason;
        return result;
      }
      stepper.update_old();
      continue;
    }

    // Linear B-field interpolation has slope jumps. Splitting at table knots
    // prevents those harmless jumps from dominating adaptive quadrature error.
    double muz_lo = sample[0].muz, muz_hi = sample[0].muz;
    for (auto const &g : sample) {
      muz_lo = std::min(muz_lo, g.muz);
      muz_hi = std::max(muz_hi, g.muz);
    }
    // muz(psi)=e1.z*cos(psi)+e2.z*sin(psi). Include its exact angular extrema,
    // even when they lie between probes, before selecting field-table knots.
    double const psi_lo = state(0)[1], psi_hi = state(1)[1];
    double const phase = std::atan2(photon.e2.z, photon.e1.z);
    for (int k = int(std::ceil((psi_lo - phase) / pi)); k <= int(std::floor((psi_hi - phase) / pi));
         ++k) {
      double const target = phase + k * pi;
      // ceil/floor of (psi-phase)/pi may round onto an adjacent integer
      // when psi is tiny (e.g. sin(pi) on a radial inward ray). Verify
      // the actual bracket before treating that integer as an extremum.
      if (options.guard_extremum_bounds && (target < psi_lo || target > psi_hi)) continue;
      auto psi_event = [&](double s) { return state(s)[1] - target; };
      double const s = locate(psi_event, 0, 1);
      auto const g = geometry(s);
      muz_lo = std::min(muz_lo, g.muz);
      muz_hi = std::max(muz_hi, g.muz);
      cuts.push_back({s});
    }
    double const dmu = (photon.bfield.mu_max - photon.bfield.mu_min) / (photon.bfield.mu_num - 1);
    if (options.fast_split_field_knots && dmu > 0 && muz_hi > muz_lo) {
      for (double sign : {-1., 1.}) {
        double const lo = sign > 0 ? muz_lo : -muz_hi;
        double const hi = sign > 0 ? muz_hi : -muz_lo;
        int const first = std::max(0, int(std::ceil((lo - photon.bfield.mu_min) / dmu)));
        int const last =
            std::min(photon.bfield.mu_num - 1, int(std::floor((hi - photon.bfield.mu_min) / dmu)));
        for (int i = first; i <= last; ++i) {
          double const knot = sign * (photon.bfield.mu_min + i * dmu);
          find_surface({Surface::magnetic_knot, knot});
        }
      }
    }
    std::sort(cuts.begin(), cuts.end(), [](Cut a, Cut b) { return a.s < b.s; });
    std::vector<Cut> unique;
    for (auto cut : cuts) {
      if (cut.s > stop) continue;
      if (unique.empty() || cut.s - unique.back().s >= 32 * std::numeric_limits<double>::epsilon())
        unique.push_back(cut);
    }
    if (unique.back().s != stop) unique.push_back({stop});

    for (std::size_t j = 1; j < unique.size(); ++j) {
      double const a = unique[j - 1].s, b = unique[j].s;
      if (!(b > a)) continue;
      if (!in_population(geometry(.5 * (a + b)))) continue;
      // Direct adaptive integration in the dense-step coordinate. No
      // singularity-specific substitution or fitted discriminant is used.
      auto segment_integral = [&](double upper) {
        if (!(upper > a)) return quad::Value{0, 0};
        double const atol = options.quadrature_atol * h * (upper - a) / (4 * r_escape);
        quad::Statistics work;
        auto value = quad::integrate(
            [&](double s) { return h * rate(geometry(s)); }, a, upper,
            {atol, options.quadrature_rtol, options.max_quadrature_intervals}, work);
        result.stats.quadrature_intervals += work.intervals;
        return value;
      };
      auto const integral = segment_integral(b);
      if (result.tau + integral.integral >= tau_target) {
        double const wanted = tau_target - result.tau;
        double lo = a, hi = b;
        // Safeguarded Newton inversion converges in a handful of integrals;
        // stopping in optical depth matches the accuracy of the integral.
        double s = a + (b - a) * wanted / integral.integral;
        bool converged = false;
        for (int it = 0; it < 48; ++it) {
          double residual = segment_integral(s).integral - wanted;
          if (std::abs(residual) <= options.quadrature_atol + options.quadrature_rtol * wanted) {
            converged = true;
            break;
          }
          if (residual < 0)
            lo = s;
          else
            hi = s;
          double derivative = h * rate(geometry(s));
          double next = s - residual / derivative;
          if (!(next > lo && next < hi) || !std::isfinite(next)) next = .5 * (lo + hi);
          if (next == s) break;
          s = next;
        }
        if (!converged) throw std::runtime_error("fast scattering inversion did not converge");
        result.distance = stepper.x_old + s * h;
        result.state = state(s);
        result.tau = result.state[3] = tau_target;
        result.termination = TransportTermination::scattered;
        return result;
      }
      result.tau += integral.integral;
      result.stats.estimated_quadrature_error += integral.error;
    }
    result.distance = stepper.x_old + stop * h;
    result.state = state(stop);
    result.state[3] = result.tau;
    if (stopped) {
      result.termination = stop_reason;
      return result;
    }
    stepper.update_old();
  }
  throw std::runtime_error("resonance transport exceeded geometry step limit");
}

// Adapter for the current repository's Photon representation.
inline TransportResult propagate(BField const &field, Boltzmann const &fb, Photon const &p,
                                 TransportOptions const &options = {},
                                 double target = std::numeric_limits<double>::infinity(),
                                 double initial_tau = 0, double escape_radius = 10000) {
  YVector state{p.r, p.psi, p.alpha, initial_tau};
  PhotonEvolution physics{field, fb, p.n, p.e1, p.e2, state, p.omega_inf, p.pol};
  return transport_fast(physics, state, options, target, escape_radius);
}
} // namespace fast_transport
