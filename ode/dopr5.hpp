#pragma once

#include "roots.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

template <int N> struct StepperDopr5 {
  static constexpr double EPS = std::numeric_limits<double>::epsilon();
  using YVector = std::array<double, N>;
  using DerivFunc = void (*)(double, YVector const &, YVector &);
  using EventFunc = double (*)(double, YVector const &);

  DerivFunc derivs;
  double x_old, h_old, h_new;
  YVector y_old, y_new, y_err;
  YVector dydx_old, dydx_new;
  double atol, rtol;

  YVector k2, k3, k4, k5, k6;
  YVector rcont1, rcont2, rcont3, rcont4, rcont5;

  double errold;
  bool reject;

  int event_count;
  std::vector<EventFunc> event_funcs;
  std::vector<int> event_signs;

  StepperDopr5(DerivFunc derivs, double atol, double rtol)
      : derivs(derivs), atol(atol), rtol(rtol) {}

  void add_event(EventFunc event_func) {
    ++event_count;
    event_funcs.push_back(event_func);
    event_signs.push_back(0);
  }

  void init(double x_init, double h_init, YVector const &y_init) {
    x_old = x_init;
    h_old = h_init;
    y_old = y_init;
    derivs(x_old, y_old, dydx_old);
    errold = 1.0e-4;
    reject = false;
    for (int i = 0; i < event_count; ++i) {
      double event_value = event_funcs[i](x_old, y_old);
      event_signs[i] = (event_value > 0) - (event_value < 0);
    }
  }

  void try_step() {
    constexpr double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, a21 = 0.2, a31 = 3.0 / 40.0,
                     a32 = 9.0 / 40.0, a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0,
                     a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                     a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                     a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
                     a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                     a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0, e1 = 71.0 / 57600.0,
                     e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0,
                     e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
    YVector y_tmp;
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * a21 * dydx_old[i];
    derivs(x_old + c2 * h_old, y_tmp, k2);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a31 * dydx_old[i] + a32 * k2[i]);
    derivs(x_old + c3 * h_old, y_tmp, k3);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a41 * dydx_old[i] + a42 * k2[i] + a43 * k3[i]);
    derivs(x_old + c4 * h_old, y_tmp, k4);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a51 * dydx_old[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    derivs(x_old + c5 * h_old, y_tmp, k5);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a61 * dydx_old[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
                                     a65 * k5[i]);
    double x_new = x_old + h_old;
    derivs(x_new, y_tmp, k6);
    for (int i = 0; i < N; i++)
      y_new[i] = y_old[i] + h_old * (a71 * dydx_old[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] +
                                     a76 * k6[i]);
    derivs(x_new, y_new, dydx_new);
    for (int i = 0; i < N; i++) {
      y_err[i] = h_old * (e1 * dydx_old[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] +
                          e7 * dydx_new[i]);
    }
  }

  void do_step() {
    while (true) {
      try_step();
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

  double error() const {
    double err = 0.0;
    for (int i = 0; i < N; i++) {
      double sk = atol + rtol * std::max(std::abs(y_old[i]), std::abs(y_new[i]));
      err += (y_err[i] / sk) * (y_err[i] / sk);
    }
    return std::sqrt(err / N);
  }

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

  void prepare_dense() {
    constexpr double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
                     d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0,
                     d6 = -1453857185.0 / 822651844.0, d7 = 69997945.0 / 29380423.0;
    for (int i = 0; i < N; i++) {
      rcont1[i] = y_old[i];
      double ydiff = y_new[i] - y_old[i];
      rcont2[i] = ydiff;
      double bspl = h_old * dydx_old[i] - ydiff;
      rcont3[i] = bspl;
      rcont4[i] = ydiff - h_old * dydx_new[i] - bspl;
      rcont5[i] = h_old * (d1 * dydx_old[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] +
                           d7 * dydx_new[i]);
    }
  }

  YVector dense_out(double x) const {
    YVector ret;
    double s = (x - x_old) / h_old;
    double s1 = 1.0 - s;
    for (int i = 0; i < N; i++)
      ret[i] = rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * rcont5[i])));
    return ret;
  }

  bool detect_event() {
    // return true if ODE system should terminate due to an event occur
    std::vector<std::pair<int, double>> active_events;
    active_events.reserve(event_count);
    double x_new = x_old + h_old;
    for (int i = 0; i < event_count; ++i) {
      double event_value = event_funcs[i](x_new, y_new);
      int event_sign_new = (event_value > 0) - (event_value < 0);
      if (event_sign_new * event_signs[i] < 0) {
        active_events.push_back(std::make_pair(i, 0.0));
      }
    }
    if (!active_events.empty()) {
      prepare_dense();
      for (auto &event : active_events) {
        event.second = zriddr([&](double x) { return event_funcs[event.first](x, dense_out(x)); },
                              x_old, x_new, 4 * EPS * (std::abs(x_old) + 1));
      }
      std::sort(active_events.begin(), active_events.end(),
                [](std::pair<int, double> const &a, std::pair<int, double> const &b) {
                  return a.second < b.second;
                });
      int event_id = active_events.front().first;
      double event_x = active_events.front().second;
      h_old = event_x - x_old;
      y_new = dense_out(event_x);
      derivs(event_x, y_new, dydx_new);

      for (int i = 0; i < event_count; ++i) {
        if (i == event_id) {
          event_signs[i] = 0;
        } else {
          double event_value = event_funcs[i](event_x, y_new);
          event_signs[i] = (event_value > 0) - (event_value < 0);
        }
      }

      return true;
    }
    return false;
  }
};
