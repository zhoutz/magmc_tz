#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

template <int N, class D> struct StepperDopr5 {
  using YVector = std::array<double, N>;
  static constexpr double EPS = std::numeric_limits<double>::epsilon();

  D &derivs;
  double x_old, h_old, h_new;
  YVector y_old, y_new, y_err;
  YVector dydx_old, dydx_new;
  double atol, rtol;

  StepperDopr5(D &derivs, double x_init, double h_init, YVector y_init, YVector dydx_init,
               double atol, double rtol)
      : derivs(derivs), x_old(x_init), h_old(h_init), y_old(y_init), dydx_old(dydx_init),
        atol(atol), rtol(rtol) {}

  static void try_step(D &derivs, double x_old, double h, YVector const &y_old,
                       YVector const &dydx_old, YVector &y_new, YVector &dydx_new, YVector &y_err) {
    constexpr double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, a21 = 0.2, a31 = 3.0 / 40.0,
                     a32 = 9.0 / 40.0, a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0,
                     a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                     a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                     a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
                     a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                     a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0, e1 = 71.0 / 57600.0,
                     e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0,
                     e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
    YVector y_tmp, k2, k3, k4, k5, k6;
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h * a21 * dydx_old[i];
    derivs(x_old + c2 * h, y_tmp, k2);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h * (a31 * dydx_old[i] + a32 * k2[i]);
    derivs(x_old + c3 * h, y_tmp, k3);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h * (a41 * dydx_old[i] + a42 * k2[i] + a43 * k3[i]);
    derivs(x_old + c4 * h, y_tmp, k4);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h * (a51 * dydx_old[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    derivs(x_old + c5 * h, y_tmp, k5);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] +
                 h * (a61 * dydx_old[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    double x_new = x_old + h;
    derivs(x_new, y_tmp, k6);
    for (int i = 0; i < N; i++)
      y_new[i] = y_old[i] +
                 h * (a71 * dydx_old[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    derivs(x_new, y_new, dydx_new);
    for (int i = 0; i < N; i++) {
      y_err[i] = h * (e1 * dydx_old[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] +
                      e7 * dydx_new[i]);
    }
  }

  void do_step() {
    while (true) {
      try_step(derivs, x_old, h_old, y_old, dydx_old, y_new, dydx_new, y_err);
      if (success(error())) break;
      if (std::abs(h_old) <= std::abs(x_old) * EPS)
        throw std::runtime_error("stepsize underflow in StepperDopr5");
    }
  }

  void update_old() {
    dydx_old = dydx_new;
    y_old = y_new;
    x_old += h_old;
    h_old = h_new;
  }

  double error() {
    double err = 0.0;
    for (int i = 0; i < N; i++) {
      double sk = atol + rtol * std::max(std::abs(y_old[i]), std::abs(y_new[i]));
      err += (y_err[i] / sk) * (y_err[i] / sk);
    }
    return std::sqrt(err / N);
  }

  double errold = 1.0e-4;
  bool reject = false;

  bool success(double err) {
    constexpr double beta = 0.4 / 5.0, alpha = 0.2 - beta * 0.75, safe = 0.9, minscale = 0.2,
                     maxscale = 10.0;
    double scale;
    if (err <= 1.0) {
      if (err == 0.0)
        scale = maxscale;
      else {
        scale = safe * std::pow(err, -alpha) * std::pow(errold, beta);
        if (scale < minscale) scale = minscale;
        if (scale > maxscale) scale = maxscale;
      }
      if (reject)
        h_new = h_old * std::min(scale, 1.0);
      else
        h_new = h_old * scale;
      errold = std::max(err, 1.0e-4);
      reject = false;
      return true;
    } else {
      scale = std::max(safe * std::pow(err, -alpha), minscale);
      h_old *= scale;
      reject = true;
      return false;
    }
  }
};
