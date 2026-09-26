// From 4be9c249402dab313c171bc04244bc2be6a39563.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>

// Adaptive open-node Gauss--Kronrod (7,15), shared by all scalar quadratures.
// Finite intervals, signed integrands and reversed bounds are supported.
// The error estimate concerns sampled function values; it cannot detect an
// entirely missed feature or repair cancellation inside the integrand.
namespace quad {
struct Options {
  double atol = 1e-9;
  double rtol = 1e-7;
  std::size_t max_intervals = 8192; // maximum active intervals
};
struct Value {
  double integral = 0, error = 0;
};
struct Statistics {
  std::size_t intervals = 0; // total GK panels evaluated, including parents
  std::size_t evaluations = 0;
};
namespace detail {
template <class Function>
Value kronrod15(Function const &fn, double a, double b, Statistics &stats) {
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
  ++stats.intervals;
  auto evaluate = [&](double x) {
    ++stats.evaluations;
    return fn(x);
  };
  double const mid = .5 * (a + b), half = .5 * (b - a);
  double const center = evaluate(mid);
  double qg = wg[3] * center, qk = wk[7] * center, qabs = wk[7] * std::abs(center);
  std::array<double, 7> fl{}, fr{};
  for (int j = 0; j < 7; ++j) {
    double const left = mid - half * x[j], right = mid + half * x[j];
    fl[j] = evaluate(left);
    fr[j] = evaluate(right);
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
  if (!std::isfinite(qk) || !std::isfinite(error))
    throw std::runtime_error("non-finite quadrature");
  return {qk, error};
}

} // namespace detail

// Statistics are accumulated into the caller's object, including failed calls.
// Accept when estimated_error <= atol + rtol * abs(integral).
template <class Function>
Value integrate(Function const &fn, double a, double b, Options const &options, Statistics &stats) {
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(options.atol) ||
      !std::isfinite(options.rtol) || options.atol < 0 || options.rtol < 0 ||
      (options.atol == 0 && options.rtol == 0) || options.max_intervals == 0)
    throw std::invalid_argument("invalid adaptive quadrature bounds or options");
  if (a == b) return {};
  bool const reverse = a > b;
  if (reverse) std::swap(a, b);
  struct Interval {
    double a, b;
    Value q;
    bool operator<(Interval const &other) const { return q.error < other.q.error; }
  };
  auto initial = detail::kronrod15(fn, a, b, stats);
  if (initial.error <= options.atol + options.rtol * std::abs(initial.integral)) {
    if (reverse) initial.integral = -initial.integral;
    return initial;
  }
  std::priority_queue<Interval> intervals;
  intervals.push({a, b, initial});
  double sum = initial.integral, error = initial.error;
  std::size_t count = 1;
  while (error > options.atol + options.rtol * std::abs(sum)) {
    if (++count > options.max_intervals)
      throw std::runtime_error("adaptive quadrature exceeded interval limit");
    Interval const top = intervals.top();
    intervals.pop();
    double const m = .5 * (top.a + top.b);
    if (m == top.a || m == top.b)
      throw std::runtime_error("adaptive quadrature interval underflow");
    auto left = detail::kronrod15(fn, top.a, m, stats);
    auto right = detail::kronrod15(fn, m, top.b, stats);
    sum += left.integral + right.integral - top.q.integral;
    error = std::max(0., error + left.error + right.error - top.q.error);
    intervals.push({top.a, m, left});
    intervals.push({m, top.b, right});
  }
  return {reverse ? -sum : sum, error};
}

template <class Function>
Value integrate(Function const &fn, double a, double b, Options const &options = {}) {
  Statistics stats;
  return integrate(fn, a, b, options, stats);
}

// Convenience overload for the transport panels and distribution quantiles.
template <class Function>
double adaptive(Function const &fn, double a, double b, double atol, double rtol,
                long &evaluations) {
  Statistics stats;
  try {
    auto value = integrate(fn, a, b, Options{atol, rtol}, stats);
    evaluations += static_cast<long>(stats.evaluations);
    return value.integral;
  } catch (...) {
    evaluations += static_cast<long>(stats.evaluations);
    throw;
  }
}
} // namespace quad
