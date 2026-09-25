#pragma once

#if 0

#include <cmath>
#include <stdexcept>

template <class F> double zriddr(F const &func, double x1, double x2, double xacc) {
  constexpr int MAXIT = 60;
  double fl = func(x1);
  double fh = func(x2);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    double xl = x1;
    double xh = x2;
    double ans = -9.99e99;
    for (int j = 0; j < MAXIT; j++) {
      double xm = 0.5 * (xl + xh);
      double fm = func(xm);
      double s = std::sqrt(fm * fm - fl * fh);
      if (s == 0.0) return ans;
      double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (std::abs(xnew - ans) <= xacc) return ans;
      ans = xnew;
      double fnew = func(ans);
      if (fnew == 0.0) return ans;
      if (std::copysign(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (std::copysign(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (std::copysign(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else {
        throw std::runtime_error("never get here.");
      }
      if (std::abs(xh - xl) <= xacc) return ans;
    }
    throw std::runtime_error("zriddr exceed maximum iterations");
  } else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    throw std::runtime_error("root must be bracketed in zriddr.");
  }
}

#else

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

template <class F>
double zriddr(F const &func, double x1, double x2, double xacc) {
  constexpr int MAXIT = 128;

  double xl = std::min(x1, x2);
  double xh = std::max(x1, x2);
  double fl = func(xl);
  double fh = func(xh);

  if (fl == 0.0)
    return xl;
  if (fh == 0.0)
    return xh;

  if (std::signbit(fl) == std::signbit(fh))
    throw std::runtime_error("Root must be bracketed");

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
    if (xm == xl || xm == xh || std::max(xm - xl, xh - xm) <= xacc)
      return xm;

    const double fm = func(xm);
    if (fm == 0.0)
      return xm;

    const double scale = std::max({std::abs(fl), std::abs(fh), std::abs(fm)});

    const double a = fl / scale;
    const double b = fh / scale;
    const double c = fm / scale;

    // fl、fh 异号，所以 sqrt(c*c - a*b)
    // 等于 hypot(c, sqrt(|a|)*sqrt(|b|))。
    const double s =
        std::hypot(c, std::sqrt(std::abs(a)) * std::sqrt(std::abs(b)));

    double xn = xm;
    if (s > 0.0) {
      xn = xm + (xm - xl) * std::copysign(1.0, fl) * (c / s);
    }

    // 每轮至少完成一次二分缩区间。
    update_bracket(xm, fm);

    // 候选点严格位于剩余括根区间内才使用。
    if (xl < xn && xn < xh) {
      const double fn = func(xn);
      if (fn == 0.0)
        return xn;

      update_bracket(xn, fn);
    }
  }

  throw std::runtime_error("zriddr exceeded maximum iterations");
}

#endif

template <class F, class DF>
double rtsafe(F const &func, DF const &dfunc, double x1, double x2,
              double xacc) {
  constexpr int MAXIT = 100;
  double xh, xl;
  double fl = func(x1);
  double fh = func(x2);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    throw std::runtime_error("Root must be bracketed in rtsafe");
  }
  if (fl == 0.0)
    return x1;
  if (fh == 0.0)
    return x2;
  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    xh = x1;
    xl = x2;
  }
  double rts = std::midpoint(x1, x2);
  double dxold = std::abs(x2 - x1);
  double dx = dxold;
  double f = func(rts);
  double df = dfunc(rts);
  for (int j = 0; j < MAXIT; j++) {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
        (std::abs(2.0 * f) > std::abs(dxold * df))) {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    } else {
      dxold = dx;
      dx = f / df;
      double temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (std::abs(dx) < xacc)
      return rts;
    f = func(rts);
    df = dfunc(rts);
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  throw std::runtime_error("Maximum number of iterations exceeded in rtsafe");
}
