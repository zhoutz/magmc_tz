#include <algorithm>
#include <cmath>
#include <limits>

inline bool solve_quadratic(double a, double b, double c, double &x1, double &x2) noexcept {
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();

  auto fail = [&]() noexcept {
    x1 = nan;
    x2 = nan;
    return false;
  };

  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)) {
    return fail();
  }

  if (a == 0.0) {
    return fail();
  }

  const double scale = std::max({std::abs(a), std::abs(b), std::abs(c)});
  if (scale == 0.0 || !std::isfinite(scale)) {
    return fail();
  }

  double A = a / scale;
  double B = b / scale;
  double C = c / scale;

  const double D = std::fma(-4.0 * A, C, B * B);

  if (!(D > 0.0)) {
    return fail();
  }

  const double sqrtD = std::sqrt(D);
  const double q = -0.5 * (B + std::copysign(sqrtD, B));

  if (q == 0.0 || !std::isfinite(q)) {
    return fail();
  }

  double r1 = q / A;
  double r2 = C / q;

  if (!std::isfinite(r1) || !std::isfinite(r2)) {
    return fail();
  }

  if (r2 < r1) {
    std::swap(r1, r2);
  }

  x1 = r1;
  x2 = r2;
  return true;
}
