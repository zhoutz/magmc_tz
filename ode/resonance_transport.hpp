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
#include <type_traits>
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
  // Fast mode omits only the exp(-u^2) Boltzmann tail beyond this landmark.
  // +infinity disables this optimization for convergence checks.
  double fast_tail_u = 12;
  bool fast_split_field_knots = false;
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

template<bool Fast = false, class Function>
QuadratureValue integrate(Function const &fn, double atol, TransportOptions const &options,
                          TransportStats &stats) {
  struct Interval {
    double a, b;
    QuadratureValue q;
    bool operator<(Interval const &other) const { return q.error < other.q.error; }
  };
  auto initial = kronrod15(fn, 0, 1, stats);
  if constexpr (Fast)
    if (initial.error <= atol + options.quadrature_rtol * std::abs(initial.integral))
      return initial;
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
template<bool Fast>
inline TransportResult transport_impl(PhotonEvolution const &photon, YVector initial,
                                 TransportOptions const &options = {},
                                 double tau_target = std::numeric_limits<double>::infinity(),
                                 double r_escape = 1000 * R_star) {
  using namespace resonance_transport_detail;
  if constexpr (Fast) if (!(options.fast_tail_u >= 8))
    throw std::invalid_argument("fast_tail_u must be at least 8 (or +infinity)");
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
  ++result.stats.field_evaluations;
  auto const radial_geometry = photon.geometry(initial, Fast);
  // Restart at the local geometry scale. Resonance splitting already protects
  // the optical-depth layer; two tiny startup steps after every scatter add no
  // accuracy. The geometry RK error controller can reject this first proposal.
  double first_step = Fast ? options.max_step_fraction * initial[0] :
                                   std::min(.01 * R_star, options.max_step_fraction * initial[0]);
  if constexpr (Fast) if (std::isfinite(tau_target)) {
    double const initial_rate = photon.rate(radial_geometry);
    if (initial_rate > 0)
      first_step = std::min(first_step, 1.5*(tau_target-initial[3])/initial_rate);
  }
  stepper.init(0, first_step,
               {initial[0], initial[1], initial[2]});

  // The reference path uses landmarks without clipping. Fast mode adds an
  // explicit, configurable Boltzmann-tail boundary for skipping empty regions.
  std::vector<Surface> surfaces{{Surface::discriminant, 0}, {Surface::beta, 0},
                                {Surface::radius, R_star}, {Surface::radius, r_escape}};
  for (double t : {.25, .5, 1., 2., 4., 8.}) {
    if constexpr (Fast) {
      if (t == .25 || t == .5 || t == 2) continue;
      if (t == 8 && std::isfinite(options.fast_tail_u)) t = options.fast_tail_u;
    }
    double const gm1 = t * t / photon.fb.a;
    double const beta = std::copysign(std::sqrt(gm1 * (gm1 + 2)) / (1 + gm1), photon.fb.b0);
    if (std::abs(beta) < 1) surfaces.push_back({Surface::beta, beta});
  }

  // A radial ray keeps its angular field and direction fixed: do not rebuild
  // a spherical basis or interpolate the same field at every quadrature node.
  bool const radial = initial[2] == 0 || initial[2] == pi;
  double const tail_gamma = std::isfinite(options.fast_tail_u) ?
                            1 + options.fast_tail_u*options.fast_tail_u/photon.fb.a : 1;
  double const tail_beta = std::copysign(std::sqrt((tail_gamma-1)*(tail_gamma+1))/tail_gamma,
                                       photon.fb.b0);
  auto in_population = [&](ResonanceGeometry const &g) {
    if (!(g.discriminant > 0)) return false;
    if (!std::isfinite(options.fast_tail_u)) return true;
    double const edge = tail_gamma*(1-tail_beta*g.mu);
    double lower = std::min(1., edge), upper = std::max(1., edge);
    if (tail_beta*g.mu > 0 && std::abs(g.mu) < std::abs(tail_beta))
      lower = std::sqrt((1-g.mu)*(1+g.mu));
    return g.x >= lower && g.x <= upper;
  };
  auto locate = [&](auto const &fn, double a, double b, double fa, double fb) {
    if constexpr (Fast) return zriddr(fn, a, b, 2 * std::numeric_limits<double>::epsilon());
    else return root(fn, a, b, fa, fb);
  };
  auto rate = [&](ResonanceGeometry const &g) {
    ++result.stats.rate_evaluations;
    if constexpr (Fast) {
      std::array<double, 2> betas{};
      auto weights = photon.resonance_weights(g, photon.pol, betas);
      return (weights[0] + weights[1]) * (photon.bfield.p + 1) * pi * g.current /
             (std::abs(photon.fb.b_bar()) * g.r);
    }
    return photon.rate(g);
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
        y[j] = stepper.rcont1[j] + s * (stepper.rcont2[j] + s1 *
            (stepper.rcont3[j] + s * (stepper.rcont4[j] + s1 * stepper.rcont5[j])));
      return {y[0], y[1], y[2], result.tau};
    };
    auto geometry = [&](double s) {
      ++result.stats.geometry_evaluations;
      if constexpr (Fast) if (radial) {
        auto g = radial_geometry;
        g.r = state(s)[0];
        g.x *= std::pow(initial[0]/g.r, 2+photon.bfield.p) *
               std::sqrt((1-rs/g.r)/(1-rs/initial[0]));
        g.discriminant = std::fma(g.x, g.x, (g.mu-1)*(g.mu+1));
        return g;
      }
      ++result.stats.field_evaluations;
      return photon.geometry(state(s), Fast);
    };
    constexpr int probes = Fast ? 2 : 8;
    std::array<ResonanceGeometry, probes + 1> sample{}, minus{}, plus{};
    constexpr double derivative_delta = 1e-5;
    for (int i = 0; i <= probes; ++i) {
      double const s = double(i) / probes;
      sample[i] = geometry(s);
      minus[i] = Fast && i == 0 ? sample[i] : geometry(std::max(0., s - derivative_delta));
      plus[i] = Fast && i == probes ? sample[i] : geometry(std::min(1., s + derivative_delta));
    }
    std::vector<Cut> cuts{{0, false}, {1, false}};
    double stop = 1;
    bool stopped = false;
    TransportTermination stop_reason = TransportTermination::escaped;
    std::vector<double> field_knots;
    auto add_cut = [&](double s, Surface const &surface) {
      ++result.stats.event_roots;
      if (surface.kind == Surface::magnetic_knot) {
        field_knots.push_back(s);
        if constexpr (Fast) if (!options.fast_split_field_knots) return;
      }
      cuts.push_back({s, surface.kind == Surface::discriminant, surface.kind == Surface::magnetic_knot});
      if (surface.kind == Surface::radius && s > 16 * std::numeric_limits<double>::epsilon() &&
          s <= stop) {
        stop = s;
        stopped = true;
        stop_reason = surface.value == R_star ? TransportTermination::absorbed
                                              : TransportTermination::escaped;
      }
    };
    bool possible = std::any_of(sample.begin(), sample.end(),
                                 [](auto const &g) { return g.discriminant > 0; });
    auto find_surface = [&](Surface const &surface) {
      double const beta_sqrt = surface.kind == Surface::beta ?
                              std::sqrt((1-surface.value)*(1+surface.value)) : 0;
      auto value = [&](ResonanceGeometry const &g) {
        if constexpr (Fast) {
          if (surface.kind == Surface::discriminant && g.discriminant > 0) possible = true;
          if (surface.kind == Surface::beta)
            return std::fma(g.x, beta_sqrt, surface.value * g.mu - 1);
        }
        return surface(g);
      };
      auto fn = [&](double s) {
        if constexpr (Fast) {
          if (surface.kind == Surface::radius) return state(s)[0] - surface.value;
          if (surface.kind == Surface::magnetic_knot) {
            double psi = state(s)[1];
            return photon.e1.z*std::cos(psi) + photon.e2.z*std::sin(psi) - surface.value;
          }
        }
        return value(geometry(s));
      };
      for (int i = 0; i < probes; ++i) {
        double const a = double(i) / probes, b = double(i + 1) / probes;
        double const fa = value(sample[i]), fb = value(sample[i + 1]);
        double const da = value(plus[i]) - value(minus[i]);
        double const db = value(plus[i + 1]) - value(minus[i + 1]);
        std::conditional_t<Fast, std::array<double,3>, std::vector<double>> local{a, b};
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
              hi = d; d = c; fd = fc;
              c = hi - ratio * (hi - lo); fc = fn(c);
            } else {
              lo = c; c = d; fc = fd;
              d = lo + ratio * (hi - lo); fd = fn(d);
            }
          }
          double const extreme = .5 * (lo + hi);
          local = {a, extreme, b};
          local_size = 3;
          // Even without a zero this point may be a narrow near-tangent peak.
          cuts.push_back({extreme, false});
        }
        for (std::size_t j = 1; j < local_size; ++j) {
          double const l = local[j - 1], r = local[j];
          double const fl = l == a ? fa : fn(l), fr = r == b ? fb : fn(r);
          if (fl == 0) add_cut(l, surface);
          if (fr == 0) add_cut(r, surface);
          if (std::signbit(fl) != std::signbit(fr))
            add_cut(locate(fn, l, r, fl, fr), surface);
        }
      }
    };
    if constexpr (Fast) {
      find_surface(surfaces[0]); // D, including interior extrema
      find_surface(surfaces[2]); // stellar surface
      find_surface(surfaces[3]); // escape surface
      if (!possible) {
        result.distance = stepper.x_old + stop*h;
        result.state = state(stop);
        result.state[3] = result.tau;
        if (stopped) { result.termination = stop_reason; return result; }
        stepper.update_old();
        continue;
      }
      for (auto const &surface : surfaces)
        if (surface.kind == Surface::beta) find_surface(surface);
    } else for (auto const &surface : surfaces) find_surface(surface);

    if constexpr (Fast) {
      std::sort(cuts.begin(), cuts.end(), [](Cut a, Cut b) { return a.s < b.s; });
      bool active = false;
      for (std::size_t j=1; j<cuts.size(); ++j) {
        double a = cuts[j-1].s, b = std::min(stop, cuts[j].s);
        if (b>a && in_population(geometry(.5*(a+b)))) { active = true; break; }
      }
      if (!active) {
        result.distance = stepper.x_old + stop*h;
        result.state = state(stop);
        result.state[3] = result.tau;
        if (stopped) { result.termination = stop_reason; return result; }
        stepper.update_old();
        continue;
      }
    }

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
      double const s = locate(psi_event, 0, 1, psi_lo - target, psi_hi - target);
      auto const g = geometry(s);
      muz_lo = std::min(muz_lo, g.muz);
      muz_hi = std::max(muz_hi, g.muz);
      cuts.push_back({s, false, false});
    }
    double const dmu = (photon.bfield.mu_max - photon.bfield.mu_min) /
                       (photon.bfield.mu_num - 1);
    std::vector<double> fold_latitudes;
    if constexpr (Fast) if (!options.fast_split_field_knots)
      for (auto const &cut : cuts)
        if (cut.discriminant_root) fold_latitudes.push_back(geometry(cut.s).muz);
    bool const need_knots = !Fast || options.fast_split_field_knots || !fold_latitudes.empty();
    if (need_knots && dmu > 0 && muz_hi > muz_lo) {
      for (double sign : {-1., 1.}) {
        double const lo = sign > 0 ? muz_lo : -muz_hi;
        double const hi = sign > 0 ? muz_hi : -muz_lo;
        int const first = std::max(0, int(std::ceil((lo - photon.bfield.mu_min) / dmu)));
        int const last = std::min(photon.bfield.mu_num - 1,
                                 int(std::floor((hi - photon.bfield.mu_min) / dmu)));
        for (int i = first; i <= last; ++i) {
          double const knot = sign * (photon.bfield.mu_min + i * dmu);
          // Without forced table splitting only the nearest knots on each
          // side of a D root are needed to bound its local regularization fit.
          if constexpr (Fast) if (!options.fast_split_field_knots &&
              std::none_of(fold_latitudes.begin(), fold_latitudes.end(),
                          [&](double mu) { return std::abs(mu-knot) <= 1.5*dmu; })) continue;
          find_surface({Surface::magnetic_knot, knot});
        }
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
      if constexpr (Fast) if (!in_population(geometry(.5*(a+b)))) continue;
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
        for (double knot : field_knots) {
          if ((knot-origin)*direction <= 0) continue;
          if ((knot-bound)*direction < 0) bound = knot;
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
          return rate(g) * h * width * pi * sine * cosine;
        };
        double const atol = options.quadrature_atol * h * width / (4 * r_escape);
        return integrate<Fast>(integrand, atol, options, result.stats);
      };
      auto const integral = segment_integral(b);
      if (result.tau + integral.integral >= tau_target) {
        double const wanted = tau_target - result.tau;
        double lo = a, hi = b;
        if constexpr (Fast) {
          // Safeguarded Newton inversion converges in a handful of integrals;
          // stopping in optical depth matches the accuracy of the integral.
          double s = a + (b-a) * wanted / integral.integral;
          bool converged = false;
          for (int it = 0; it < 48; ++it) {
            double residual = segment_integral(s).integral - wanted;
            if (std::abs(residual) <= options.quadrature_atol +
                options.quadrature_rtol * wanted) { converged = true; break; }
            if (residual < 0) lo = s; else hi = s;
            double derivative = h * rate(geometry(s));
            double next = s - residual/derivative;
            if (!(next > lo && next < hi) || !std::isfinite(next)) next = .5*(lo+hi);
            if (next == s) break;
            s = next;
          }
          if (!converged) throw std::runtime_error("fast scattering inversion did not converge");
          result.distance = stepper.x_old + s*h;
          result.state = state(s);
          result.tau = result.state[3] = tau_target;
          result.termination = TransportTermination::scattered;
          return result;
        }
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

// Retain the original implementation for reproducible speed/accuracy comparisons.
inline TransportResult transport(PhotonEvolution const &photon, YVector initial,
                                 TransportOptions const &options = {},
                                 double target = std::numeric_limits<double>::infinity(),
                                 double outer = 1000*R_star) {
  return transport_impl<false>(photon, initial, options, target, outer);
}

inline TransportResult transport_fast(PhotonEvolution const &photon, YVector initial,
                                      TransportOptions const &options = {},
                                      double target = std::numeric_limits<double>::infinity(),
                                      double outer = 1000*R_star) {
  return transport_impl<true>(photon, initial, options, target, outer);
}
