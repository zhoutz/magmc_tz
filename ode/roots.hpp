#pragma once

#include <cmath>
#include <stdexcept>

template <class T> double zriddr(T &func, double x1, double x2, double xacc) {
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
