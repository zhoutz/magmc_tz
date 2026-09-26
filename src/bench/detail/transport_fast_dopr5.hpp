// From 4be9c249402dab313c171bc04244bc2be6a39563.
#pragma once

#include "../../roots.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

// Historical stepper kept isolated to reproduce the requested commit.
namespace fast_transport {
template <int N, class DerivFunc> struct StepperDopr5 {
  static_assert(N > 0, "StepperDopr5 requires at least one state variable");
  static constexpr double EPS = std::numeric_limits<double>::epsilon();
  using YVector = std::array<double, N>;
  using EventFunc = std::function<double(double, YVector const &)>;
  // Return a positive maximum absolute step; +infinity leaves it unrestricted.
  using StepLimitFunc = std::function<double(double, YVector const &)>;

  struct Statistics {
    std::size_t derivative_calls = 0;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::size_t event_evaluations = 0;
    std::size_t limiter_calls = 0;
  } stats;

  DerivFunc const &derivs;
  double x_old, h_old, h_new;
  YVector y_old, y_new, y_err;
  YVector dydx_old, dydx_new;
  double atol, rtol;

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
  StepLimitFunc step_limiter;
  bool initialized = false;
  bool dense_ready = false;
  double dense_x_old = 0.0, dense_h = 0.0;

  StepperDopr5(DerivFunc const &derivs, double atol, double rtol)
      : derivs(derivs), atol(atol), rtol(rtol) {
    if (!std::isfinite(atol) || !std::isfinite(rtol) || atol < 0.0 || rtol < 0.0 ||
        (atol == 0.0 && rtol == 0.0))
      throw std::invalid_argument("Invalid tolerances in StepperDopr5");
  }
  StepperDopr5(DerivFunc const &&, double, double) = delete;

  static int sign(double value) { return (value > 0) - (value < 0); }

  void reset_statistics() { stats = {}; }

  void set_step_limiter(StepLimitFunc limiter) { step_limiter = std::move(limiter); }

  static void check_finite(YVector const &y) {
    for (double value : y)
      if (!std::isfinite(value))
        throw std::runtime_error("Non-finite state or derivative in StepperDopr5");
  }

  void evaluate_derivative(double x, YVector const &y, YVector &dy) {
    if (!std::isfinite(x)) throw std::runtime_error("Non-finite abscissa in StepperDopr5");
    check_finite(y);
    ++stats.derivative_calls;
    derivs(x, y, dy);
    check_finite(dy);
  }

  double evaluate_event(Event const &event, double x, YVector const &y) {
    ++stats.event_evaluations;
    double value = event.func(x, y);
    if (!std::isfinite(value))
      throw std::runtime_error("Non-finite event value in StepperDopr5");
    return value;
  }

  void add_event(EventFunc event_func) {
    if (!event_func) throw std::invalid_argument("Empty event callback in StepperDopr5");
    events.emplace_back(Event{
        .func = std::move(event_func),
        .x = 0,
        .sign = 0,
        .active = false,
    });
    if (initialized) events.back().sign = sign(evaluate_event(events.back(), x_old, y_old));
  }

  void init(double x_init, double h_init, YVector const &y_init) {
    if (!std::isfinite(x_init) || !std::isfinite(h_init) || h_init == 0.0)
      throw std::invalid_argument("Invalid initial abscissa or step in StepperDopr5");
    check_finite(y_init);
    initialized = false;
    dense_ready = false;
    x_old = x_init;
    h_old = h_init;
    h_new = h_init;
    y_old = y_init;
    y_new = y_init;
    evaluate_derivative(x_old, y_old, dydx_old);
    dydx_new = dydx_old;
    errold = 1.0e-4;
    reject = false;
    for (auto &event : events) {
      event.sign = sign(evaluate_event(event, x_old, y_old));
      event.active = false;
    }
    initialized = true;
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
    evaluate_derivative(x_old + c2 * h_old, y_tmp, k2);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a31 * dydx_old[i] + a32 * k2[i]);
    evaluate_derivative(x_old + c3 * h_old, y_tmp, k3);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a41 * dydx_old[i] + a42 * k2[i] + a43 * k3[i]);
    evaluate_derivative(x_old + c4 * h_old, y_tmp, k4);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a51 * dydx_old[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    evaluate_derivative(x_old + c5 * h_old, y_tmp, k5);
    for (int i = 0; i < N; i++)
      y_tmp[i] = y_old[i] + h_old * (a61 * dydx_old[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
                                     a65 * k5[i]);
    double x_new = x_old + h_old;
    evaluate_derivative(x_new, y_tmp, k6);
    for (int i = 0; i < N; i++)
      y_new[i] = y_old[i] + h_old * (a71 * dydx_old[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] +
                                     a76 * k6[i]);
    evaluate_derivative(x_new, y_new, dydx_new);
    for (int i = 0; i < N; i++) {
      y_err[i] = h_old * (e1 * dydx_old[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] +
                          e7 * dydx_new[i]);
    }
  }

  void do_step() {
    if (!initialized) throw std::logic_error("StepperDopr5 must be initialized before stepping");
    dense_ready = false;
    if (step_limiter) {
      ++stats.limiter_calls;
      const double limit = step_limiter(x_old, y_old);
      if (!(limit > 0.0))
        throw std::runtime_error("Step limiter must return a positive bound in StepperDopr5");
      h_old = std::copysign(std::min(std::abs(h_old), limit), h_old);
    }
    while (true) {
      if (!std::isfinite(h_old) || !std::isfinite(x_old + h_old))
        throw std::runtime_error("Non-finite step or endpoint in StepperDopr5");
      if (x_old + h_old == x_old) {
        throw std::runtime_error("stepsize underflow in StepperDopr5");
      }
      try_step();
      if (success(error())) {
        ++stats.accepted_steps;
        break;
      }
      ++stats.rejected_steps;
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
        throw std::runtime_error("Non-finite error encountered in StepperDopr5");
      }
      // A component fixed exactly at zero needs no absolute tolerance.
      const double scaled_error = sk == 0.0 && y_err[i] == 0.0 ? 0.0 : std::abs(y_err[i]) / sk;
      if (!std::isfinite(scaled_error))
        throw std::runtime_error("Non-finite scaled error in StepperDopr5");
      err = std::max(err, scaled_error);
    }
    return err;
  }

  bool success(double err) {
    if (!std::isfinite(err) || err < 0.0)
      throw std::runtime_error("Invalid error estimate in StepperDopr5");
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
    // Event truncation changes h_old and y_new, but the accepted RK stages still
    // describe the full original step. Retain that polynomial and its interval.
    if (dense_ready) return;
    dense_x_old = x_old;
    dense_h = h_old;
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
    dense_ready = true;
  }

  YVector dense_out(double x) const {
    if (!dense_ready) throw std::logic_error("Call prepare_dense before dense_out");
    YVector ret;
    double s = (x - dense_x_old) / dense_h;
    double s1 = 1.0 - s;
    for (int i = 0; i < N; i++)
      ret[i] = rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * rcont5[i])));
    return ret;
  }

  int detect_event() {
    double x_new = x_old + h_old;
    for (auto &event : events) {
      double event_value = evaluate_event(event, x_new, y_new);
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
              zriddr([&](double x) { return evaluate_event(event, x, dense_out(x)); }, x_old,
                     x_new, x_tol);
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
      evaluate_derivative(leftmost_x, y_new, dydx_new);
      h_old = leftmost_x - x_old;

      for (auto &event : events) {
        if (event.active && std::abs(event.x - leftmost_x) < x_tol) {
          event.sign = 0;
        } else {
          event.sign = sign(evaluate_event(event, leftmost_x, y_new));
          event.active = false;
        }
      }

      return std::distance(events.begin(), leftmost_event);
    }
    return -1;
  }
};

}
