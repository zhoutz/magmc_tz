#pragma once

// Propagate only the three geometrical variables with DOPRI5. Optical depth is
// integrated separately on dense trajectories, split at resonance surfaces.
#include "dopr5.hpp"
#include "transport_physics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

struct TransportOptions {
  double geometry_rtol = 1e-9;
  double geometry_atol = 1e-10;
  double quadrature_rtol = 1e-7;
  double quadrature_atol = 1e-9;
  double max_step_fraction = 0.1;
  // Fit scale in normalized geometrical-step coordinates; expose it so the
  // root-regularization bias can be checked independently of GK tolerance.
  double discriminant_fit_step = 2.5e-4;
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

// Open-node Gauss--Kronrod (7,15). In particular, rate is never sampled at a
// root where the E-mode rate has an integrable inverse-square-root singularity.
struct QuadratureValue { double integral, error; };
template<class Function>
QuadratureValue kronrod15(Function const &fn, double a, double b,
                         TransportStats &stats) {
  constexpr std::array<double, 8> x{
      .991455371120812639206854697526329, .949107912342758524526189684047851,
      .864864423359769072789712788640926, .741531185599394439863864773280788,
      .586087235467691130294144838258730, .405845151377397166906606412076961,
      .207784955007898467600689403773245, 0};
  constexpr std::array<double, 8> wk{
      .022935322010529224963732008058970, .063092092629978553290700663189204,
      .104790010322250183839876322541518, .140653259715525918745189590510238,
      .169004726639267902826583426598550, .190350578064785409913256402421014,
      .204432940075298892414161999234649, .209482141084727828012999174891714};
  constexpr std::array<double, 4> wg{
      .129484966168869693270611432679082, .279705391489276667901467771423780,
      .381830050505118944950369775488975, .417959183673469387755102040816327};
  ++stats.quadrature_intervals;
  double const mid = .5 * (a + b), half = .5 * (b - a);
  double const center = fn(mid);
  double qg = wg[3] * center, qk = wk[7] * center, qabs = wk[7] * std::abs(center);
  std::array<double, 7> fl{}, fr{};
  for (int j = 0; j < 7; ++j) {
    fl[j] = fn(mid - half * x[j]);
    fr[j] = fn(mid + half * x[j]);
    double const sum = fl[j] + fr[j];
    qk += wk[j] * sum;
    qabs += wk[j] * (std::abs(fl[j]) + std::abs(fr[j]));
    if (j % 2) qg += wg[(j - 1) / 2] * sum;
  }
  double const mean = .5 * qk;
  double asc = wk[7] * std::abs(center - mean);
  for (int j = 0; j < 7; ++j)
    asc += wk[j] * (std::abs(fl[j] - mean) + std::abs(fr[j] - mean));
  qk *= half;
  asc *= half;
  double error = std::abs((qk / half - qg) * half);
  if (asc > 0 && error > 0) error = asc * std::min(1., std::pow(200 * error / asc, 1.5));
  error = std::max(error, 50 * std::numeric_limits<double>::epsilon() * half * qabs);
  if (!std::isfinite(qk) || !std::isfinite(error) || qk < 0)
    throw std::runtime_error("non-finite or negative resonance quadrature");
  return {qk, error};
}

template<class Function>
QuadratureValue integrate(Function const &fn, double atol, TransportOptions const &options,
                          TransportStats &stats) {
  struct Interval {
    double a, b;
    QuadratureValue q;
    bool operator<(Interval const &other) const { return q.error < other.q.error; }
  };
  auto initial = kronrod15(fn, 0, 1, stats);
  std::priority_queue<Interval> intervals;
  intervals.push({0, 1, initial});
  double sum = initial.integral, error = initial.error;
  std::size_t count = 1;
  while (error > atol + options.quadrature_rtol * std::abs(sum)) {
    if (++count > options.max_quadrature_intervals)
      throw std::runtime_error("resonance quadrature exceeded interval limit: integral=" + std::to_string(sum) + " error=" + std::to_string(error));
    Interval const top = intervals.top();
    intervals.pop();
    double const m = .5 * (top.a + top.b);
    if (m == top.a || m == top.b)
      throw std::runtime_error("resonance quadrature interval underflow");
    auto left = kronrod15(fn, top.a, m, stats);
    auto right = kronrod15(fn, m, top.b, stats);
    sum += left.integral + right.integral - top.q.integral;
    error = std::max(0., error + left.error + right.error - top.q.error);
    intervals.push({top.a, m, left});
    intervals.push({m, top.b, right});
  }
  return {sum, error};
}

// Work in a normalized coordinate on each geometrical step; bisection then
// retains accuracy even after a long path has accumulated a large distance.
template<class Function>
double root(Function const &fn, double a, double b, double fa, double fb) {
  if (fa == 0) return a;
  if (fb == 0) return b;
  for (int i = 0; i < 64; ++i) {
    double const m = .5 * (a + b);
    if (m == a || m == b) break;
    double const fm = fn(m);
    if (fm == 0) return m;
    if (std::signbit(fa) != std::signbit(fm)) { b = m; fb = fm; }
    else { a = m; fa = fm; }
  }
  return std::abs(fa) < std::abs(fb) ? a : b;
}

struct Surface {
  enum Kind { beta, discriminant, radius, magnetic_knot } kind;
  double value;
  double operator()(ResonanceGeometry const &g) const {
    switch (kind) {
      case beta:
        return std::log(g.x) + .5 * std::log1p(-value * value) -
               std::log1p(-value * g.mu);
      case discriminant: return g.discriminant;
      case radius: return g.r - value;
      case magnetic_knot: return g.muz - value;
    }
    return 0;
  }
};

struct Cut { double s; bool discriminant_root = false; bool magnetic_knot = false; };
} // namespace resonance_transport_detail

// tau_target is an absolute accumulated optical depth; initial[3] is retained.
// No RNG is used. A caller may pass -log(U) and scatter at the returned state.
inline TransportResult transport(PhotonEvolution const &photon, YVector initial,
                                 TransportOptions const &options = {},
                                 double tau_target = std::numeric_limits<double>::infinity(),
                                 double r_escape = 1000 * R_star) {
  using namespace resonance_transport_detail;
  if (!(options.geometry_rtol > 0) || !(options.geometry_atol > 0) ||
      !(options.quadrature_rtol > 0) || !(options.quadrature_atol > 0) ||
      !(options.max_step_fraction > 0 && options.max_step_fraction <= .25) ||
      !(options.discriminant_fit_step > 0 && options.discriminant_fit_step <= .001) ||
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
  stepper.init(0, std::min(.01 * R_star, options.max_step_fraction * initial[0]),
               {initial[0], initial[1], initial[2]});

  // These are integration landmarks, not distribution cutoffs. Every interval
  // from launch to termination is still integrated, including both tails.
  std::vector<Surface> surfaces{{Surface::discriminant, 0}, {Surface::beta, 0},
                                {Surface::radius, R_star}, {Surface::radius, r_escape}};
  for (double t : {.25, .5, 1., 2., 4., 8.}) {
    double const gm1 = t * t / photon.fb.a;
    double const beta = std::copysign(std::sqrt(gm1 * (gm1 + 2)) / (1 + gm1), photon.fb.b0);
    if (std::abs(beta) < 1) surfaces.push_back({Surface::beta, beta});
  }

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
        y[j] = stepper.rcont1[j] + s * (stepper.rcont2[j] + s1 *
            (stepper.rcont3[j] + s * (stepper.rcont4[j] + s1 * stepper.rcont5[j])));
      return {y[0], y[1], y[2], result.tau};
    };
    auto geometry = [&](double s) {
      ++result.stats.geometry_evaluations;
      ++result.stats.field_evaluations;
      return photon.geometry(state(s));
    };
    constexpr int probes = 8;
    std::array<ResonanceGeometry, probes + 1> sample{}, minus{}, plus{};
    constexpr double derivative_delta = 1e-5;
    for (int i = 0; i <= probes; ++i) {
      double const s = double(i) / probes;
      sample[i] = geometry(s);
      minus[i] = geometry(std::max(0., s - derivative_delta));
      plus[i] = geometry(std::min(1., s + derivative_delta));
    }
    std::vector<Cut> cuts{{0, false}, {1, false}};
    double stop = 1;
    bool stopped = false;
    TransportTermination stop_reason = TransportTermination::escaped;
    auto add_cut = [&](double s, Surface const &surface) {
      cuts.push_back({s, surface.kind == Surface::discriminant, surface.kind == Surface::magnetic_knot});
      ++result.stats.event_roots;
      if (surface.kind == Surface::radius && s > 16 * std::numeric_limits<double>::epsilon() &&
          s <= stop) {
        stop = s;
        stopped = true;
        stop_reason = surface.value == R_star ? TransportTermination::absorbed
                                              : TransportTermination::escaped;
      }
    };
    auto find_surface = [&](Surface const &surface) {
      auto fn = [&](double s) { return surface(geometry(s)); };
      for (int i = 0; i < probes; ++i) {
        double const a = double(i) / probes, b = double(i + 1) / probes;
        double const fa = surface(sample[i]), fb = surface(sample[i + 1]);
        double const da = surface(plus[i]) - surface(minus[i]);
        double const db = surface(plus[i + 1]) - surface(minus[i + 1]);
        std::vector<double> local{a, b};
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
              hi = d; d = c; fd = fc;
              c = hi - ratio * (hi - lo); fc = fn(c);
            } else {
              lo = c; c = d; fc = fd;
              d = lo + ratio * (hi - lo); fd = fn(d);
            }
          }
          double const extreme = .5 * (lo + hi);
          local = {a, extreme, b};
          // Even without a zero this point may be a narrow near-tangent peak.
          cuts.push_back({extreme, false});
        }
        for (std::size_t j = 1; j < local.size(); ++j) {
          double const l = local[j - 1], r = local[j];
          double const fl = l == a ? fa : fn(l), fr = r == b ? fb : fn(r);
          if (fl == 0) add_cut(l, surface);
          if (fr == 0) add_cut(r, surface);
          if (std::signbit(fl) != std::signbit(fr))
            add_cut(root(fn, l, r, fl, fr), surface);
        }
      }
    };
    for (auto const &surface : surfaces) find_surface(surface);

    // Linear B-field interpolation has slope jumps. Splitting at table knots
    // prevents those harmless jumps from dominating adaptive quadrature error.
    double muz_lo = sample[0].muz, muz_hi = sample[0].muz;
    for (auto const &g : sample) { muz_lo = std::min(muz_lo, g.muz); muz_hi = std::max(muz_hi, g.muz); }
    // muz(psi)=e1.z*cos(psi)+e2.z*sin(psi). Include its exact angular extrema,
    // even when they lie between probes, before selecting field-table knots.
    double const psi_lo = state(0)[1], psi_hi = state(1)[1];
    double const phase = std::atan2(photon.e2.z, photon.e1.z);
    for (int k = int(std::ceil((psi_lo - phase) / pi));
         k <= int(std::floor((psi_hi - phase) / pi)); ++k) {
      double const target = phase + k * pi;
      auto psi_event = [&](double s) { return state(s)[1] - target; };
      double const s = root(psi_event, 0, 1, psi_lo - target, psi_hi - target);
      auto const g = geometry(s);
      muz_lo = std::min(muz_lo, g.muz);
      muz_hi = std::max(muz_hi, g.muz);
      cuts.push_back({s, false, false});
    }
    double const dmu = (photon.bfield.mu_max - photon.bfield.mu_min) /
                       (photon.bfield.mu_num - 1);
    if (dmu > 0 && muz_hi > muz_lo) {
      for (double sign : {-1., 1.}) {
        double const lo = sign > 0 ? muz_lo : -muz_hi;
        double const hi = sign > 0 ? muz_hi : -muz_lo;
        int const first = std::max(0, int(std::ceil((lo - photon.bfield.mu_min) / dmu)));
        int const last = std::min(photon.bfield.mu_num - 1,
                                 int(std::floor((hi - photon.bfield.mu_min) / dmu)));
        for (int i = first; i <= last; ++i)
          find_surface({Surface::magnetic_knot, sign * (photon.bfield.mu_min + i * dmu)});
      }
    }
    std::sort(cuts.begin(), cuts.end(), [](Cut a, Cut b) { return a.s < b.s; });
    std::vector<Cut> unique;
    for (auto cut : cuts) {
      if (cut.s > stop) continue;
      if (!unique.empty() && cut.s - unique.back().s < 32 * std::numeric_limits<double>::epsilon())
        {
        unique.back().discriminant_root |= cut.discriminant_root;
        unique.back().magnetic_knot |= cut.magnetic_knot;
      }
      else unique.push_back(cut);
    }
    if (unique.back().s != stop) unique.push_back({stop, false});

    for (std::size_t j = 1; j < unique.size(); ++j) {
      double const a = unique[j - 1].s, b = unique[j].s;
      if (!(b > a)) continue;
      auto const ga = geometry(a), gb = geometry(b);
      struct DiscriminantModel {
        bool active = false;
        double root = 0, delta = 0;
        std::array<double, 4> values{};
        double operator()(double offset) const {
          double const z = offset / delta;
          return z * (-values[0] * (z-2)*(z-3)*(z-4)/6 +
                       values[1] * (z-1)*(z-3)*(z-4)/4 -
                       values[2] * (z-1)*(z-2)*(z-4)/6 +
                       values[3] * (z-1)*(z-2)*(z-3)/24);
        }
      };
      auto make_model = [&](std::size_t index, int direction) {
        DiscriminantModel model;
        if (!unique[index].discriminant_root) return model;
        double const origin = unique[index].s;
        double bound = direction > 0 ? 1. : 0.;
        // The fit may cross thermal landmarks, but never interpolation knots.
        for (auto const &cut : unique) {
          if (!cut.magnetic_knot || (cut.s-origin)*direction <= 0) continue;
          if ((cut.s-bound)*direction < 0) bound = cut.s;
        }
        double const available = (bound-origin)*direction;
        if (!(available > 64*std::numeric_limits<double>::epsilon())) return model;
        model.active = true;
        model.root = origin;
        model.delta = direction * std::min(options.discriminant_fit_step, .25*available);
        auto const endpoint = direction > 0 ? ga : gb;
        for (int k=0; k<4; ++k) {
          auto const g = geometry(origin + (k+1)*model.delta);
          model.values[k] = (g.x-endpoint.x)*(g.x+endpoint.x) +
                            (g.mu-endpoint.mu)*(g.mu+endpoint.mu);
        }
        return model;
      };
      auto const left_model = make_model(j-1, 1), right_model = make_model(j, -1);
      auto segment_integral = [&](double upper) {
        if (!(upper > a)) return QuadratureValue{0, 0};
        double const width = upper - a;
        auto integrand = [&](double t) {
          double const angle = .5 * pi * t;
          double const sine = std::sin(angle), cosine = std::cos(angle);
          // This symmetric substitution regularizes either D=0 endpoint.
          double const left_offset = width * sine * sine;
          double const right_offset = (b - upper) + width * cosine * cosine;
          double s = t < .5 ? a + left_offset : upper - width * cosine * cosine;
          if (s <= a) s = std::nextafter(a, upper);
          if (s >= upper) s = std::nextafter(upper, a);
          auto g = geometry(s);
          // A root-anchored, one-sided quartic gives a smooth D close to
          // coalescence, where subtracting double-precision x~1 destroys
          // the sqrt(D) cancellation against the transformation Jacobian.
          if (left_model.active && left_offset < std::abs(left_model.delta) &&
              (!right_model.active || left_offset <= right_offset))
            g.discriminant = left_model(left_offset);
          else if (right_model.active && right_offset < std::abs(right_model.delta))
            g.discriminant = right_model(-right_offset);
          ++result.stats.rate_evaluations;
          return photon.rate(g) * h * width * pi * sine * cosine;
        };
        double const atol = options.quadrature_atol * h * width / (4 * r_escape);
        return integrate(integrand, atol, options, result.stats);
      };
      auto const integral = segment_integral(b);
      if (result.tau + integral.integral >= tau_target) {
        double const wanted = tau_target - result.tau;
        double lo = a, hi = b;
        // Invert the nonnegative integrated optical depth, never a DOPRI dense
        // polynomial in tau (which need not be monotone across a thin layer).
        for (int iteration = 0; iteration < 60; ++iteration) {
          double const m = .5 * (lo + hi);
          if (m == lo || m == hi) break;
          double const value = segment_integral(m).integral;
          if (value < wanted) lo = m; else hi = m;
          if (hi - lo < 8 * std::numeric_limits<double>::epsilon() * std::max(1., std::abs(m))) break;
        }
        double const s = .5 * (lo + hi);
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
    if (stopped) { result.termination = stop_reason; return result; }
    stepper.update_old();
  }
  throw std::runtime_error("resonance transport exceeded geometry step limit");
}
