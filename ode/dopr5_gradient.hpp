#pragma once

#include "roots.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

// StepperDopr5 with optical depth gradient monitoring
// Monitors dτ/dl and constrains step size to prevent missing resonances
template <int N, class DerivFunc> struct StepperDopr5Gradient {
  static_assert(N > 0, "StepperDopr5Gradient requires at least one state variable");
  static constexpr double EPS = std::numeric_limits<double>::epsilon();
  using YVector = std::array<double, N>;
  using EventFunc = std::function<double(double, YVector const &)>;

  DerivFunc const &derivs;
  double x_old, h_old, h_new;
  YVector y_old, y_new, y_err;
  YVector dydx_old, dydx_new;
  double atol, rtol;

  // Gradient monitoring parameters
  double max_dtau_per_step;  // Maximum allowed Δτ per step
  double dtau_threshold;      // Threshold for "large" dτ/dl
  int tau_index;              // Index of τ in state vector (typically 3)

  YVector k2, k3, k4, k5, k6;
  YVector rcont1, rcont2, rcont3, rcont4, rcont5;

  double errold;
  bool reject;

  struct Event {
    EventFunc func;
    double x;
    int sign;
    bool active;
  };

  std::vector<Event> events;

  StepperDopr5Gradient(DerivFunc const &derivs, double atol, double rtol,
                       double max_dtau_per_step = 0.1, double dtau_threshold = 10.0,
                       int tau_index = 3)
      : derivs(derivs), atol(atol), rtol(rtol), max_dtau_per_step(max_dtau_per_step),
        dtau_threshold(dtau_threshold), tau_index(tau_index) {}
  StepperDopr5Gradient(DerivFunc const &&, double, double, double, double, int) = delete;

  static int sign(double value) { return (value > 0) - (value < 0); }

  void add_event(EventFunc event_func) {
    events.emplace_back(Event{
        .func = std::move(event_func),
        .x = 0,
        .sign = 0,
        .active = true,
    });
  }

  void init(double x_init, double h_init, YVector const &y_init) {
    x_old = x_init;
    h_old = h_init;
    y_old = y_init;
    derivs(x_old, y_old, dydx_old);
    errold = 1.0e-4;
    reject = false;
    for (auto &event : events) {
      event.sign = sign(event.func(x_old, y_old));
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
      if (std::abs(h_old) <= std::abs(x_old) * EPS) {
        throw std::runtime_error("stepsize underflow in StepperDopr5Gradient");
      }
      try_step();
      if (success(error())) break;
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
      if (!std::isfinite(sk) || !std::isfinite(y_err[i])) {
        throw std::runtime_error("Non-finite error encountered in StepperDopr5Gradient");
      }
      err = std::max(err, std::abs(y_err[i]) / sk);
    }
    return err;
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

      // Apply gradient-based step size constraint
      double max_dtaudl = std::max({std::abs(dydx_old[tau_index]), std::abs(k2[tau_index]),
                                    std::abs(k3[tau_index]), std::abs(k4[tau_index]),
                                    std::abs(k5[tau_index]), std::abs(k6[tau_index]),
                                    std::abs(dydx_new[tau_index])});

      double h_proposed = h_old * scale;

      // If dτ/dl is large, limit step size
      if (max_dtaudl > dtau_threshold) {
        double h_gradient_limit = max_dtau_per_step / max_dtaudl;
        h_proposed = std::min(h_proposed, h_gradient_limit);
      }

      if (reject)
        h_new = std::min(h_proposed, h_old);
      else
        h_new = h_proposed;

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

  int detect_event() {
    double x_new = x_old + h_old;
    for (auto &event : events) {
      double event_value = event.func(x_new, y_new);
      int event_sign_new = sign(event_value);
      if (event.sign != 0 && event_sign_new * event.sign <= 0) {
        event.active = true;
      } else {
        event.active = false;
      }
      event.sign = event_sign_new;
    }

    if (std::any_of(events.begin(), events.end(),
                    [](Event const &event) { return event.active; })) {
      prepare_dense();
      double x_tol = 4 * EPS * (1 + std::abs(x_old));
      for (auto &event : events) {
        if (event.active) {
          event.x =
              zriddr([&](double x) { return event.func(x, dense_out(x)); }, x_old, x_new, x_tol);
        }
      }

      auto leftmost_event = events.end();
      for (auto event = events.begin(); event != events.end(); ++event) {
        if (!event->active) continue;
        if (leftmost_event == events.end() ||
            (event->x - x_old) / h_old < (leftmost_event->x - x_old) / h_old) {
          leftmost_event = event;
        }
      }
      double leftmost_x = leftmost_event->x;
      y_new = dense_out(leftmost_x);
      derivs(leftmost_x, y_new, dydx_new);
      h_old = leftmost_x - x_old;

      for (auto &event : events) {
        if (event.active && std::abs(event.x - leftmost_x) < x_tol) {
          event.sign = 0;
        } else {
          event.sign = sign(event.func(leftmost_x, y_new));
        }
      }

      return std::distance(events.begin(), leftmost_event);
    }
    return -1;
  }
};
