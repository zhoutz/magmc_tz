#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

template <class F> double zriddr(F const &func, double x1, double x2, double xacc) {
  constexpr int MAXIT = 128;

  if (!std::isfinite(x1) || !std::isfinite(x2) || !std::isfinite(xacc) || xacc < 0.0)
    throw std::invalid_argument("Invalid interval or tolerance in zriddr");
  auto evaluate = [&](double x) {
    double value = func(x);
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite function value in zriddr");
    return value;
  };

  double xl = std::min(x1, x2);
  double xh = std::max(x1, x2);
  double fl = evaluate(xl);
  double fh = evaluate(xh);

  if (fl == 0.0) return xl;
  if (fh == 0.0) return xh;

  if (std::signbit(fl) == std::signbit(fh)) throw std::runtime_error("Root must be bracketed");

  // fx 非零；保持两端函数值异号。
  auto update_bracket = [&](double x, double fx) {
    if (std::signbit(fx) == std::signbit(fl)) {
      xl = x;
      fl = fx;
    } else {
      xh = x;
      fh = fx;
    }
  };

  for (int i = 0; i < MAXIT; ++i) {
    const double xm = std::midpoint(xl, xh);

    // 区间足够小，或已无可表示的内部浮点数。
    if (xm == xl || xm == xh || std::max(xm - xl, xh - xm) <= xacc) return xm;

    const double fm = evaluate(xm);
    if (fm == 0.0) return xm;

    const double scale = std::max({std::abs(fl), std::abs(fh), std::abs(fm)});

    const double a = fl / scale;
    const double b = fh / scale;
    const double c = fm / scale;

    // fl、fh 异号，所以 sqrt(c*c - a*b)
    // 等于 hypot(c, sqrt(|a|)*sqrt(|b|))。
    const double s = std::hypot(c, std::sqrt(std::abs(a)) * std::sqrt(std::abs(b)));

    double xn = xm;
    if (s > 0.0) {
      xn = xm + (xm - xl) * std::copysign(1.0, fl) * (c / s);
    }

    // 每轮至少完成一次二分缩区间。
    update_bracket(xm, fm);

    // 候选点严格位于剩余括根区间内才使用。
    if (xl < xn && xn < xh) {
      const double fn = evaluate(xn);
      if (fn == 0.0) return xn;

      update_bracket(xn, fn);
    }
  }

  throw std::runtime_error("zriddr exceeded maximum iterations");
}
