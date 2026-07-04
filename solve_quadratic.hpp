#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

inline std::pair<double, double> solve_quadratic_real(double a, double b, double c) noexcept {
  const double nan = std::numeric_limits<double>::quiet_NaN();

  auto no_real_roots = [nan]() noexcept { return std::pair<double, double>{nan, nan}; };

  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
    return no_real_roots();
  }

  if (a == 0.0) {
    if (b == 0.0) {
      return no_real_roots();
    }

    double x = -c / b;
    return {x, x};
  }

  const double scale = std::max({std::abs(a), std::abs(b), std::abs(c)});

  using real = long double;

  const real A = static_cast<real>(a) / scale;
  const real B = static_cast<real>(b) / scale;
  const real C = static_cast<real>(c) / scale;

  const real D = std::fma(-4.0L * A, C, B * B);

  if (D < 0.0L) {
    return no_real_roots();
  }

  if (D == 0.0L) {
    double x = static_cast<double>(-B / (2.0L * A));
    return {x, x};
  }

  const real sqrtD = std::sqrt(D);
  const real q = -0.5L * (B + std::copysign(sqrtD, B));
  real x1 = q / A;
  real x2 = C / q;

  double r1 = static_cast<double>(x1);
  double r2 = static_cast<double>(x2);

  if (r2 < r1) {
    std::swap(r1, r2);
  }

  return {r1, r2};
}
