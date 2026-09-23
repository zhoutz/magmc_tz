#include "dopr5.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, std::string const &message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, std::string const &message) {
  require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
          message + ": got " + std::to_string(actual) + ", expected " + std::to_string(expected));
}

template <class F> void throws(F const &f, std::string const &message) {
  try {
    f();
  } catch (std::exception const &) {
    return;
  }
  throw std::runtime_error(message);
}

using State = std::array<double, 1>;

void test_limiter_and_statistics() {
  auto constant = [](double, State const &, State &dy) { dy[0] = 1.0; };
  StepperDopr5<1, decltype(constant)> stepper(constant, 1e-12, 1e-12);
  stepper.set_step_limiter([](double, State const &) { return 0.125; });
  stepper.init(0.0, 1.0, {0.0});
  stepper.do_step();
  near(stepper.h_old, 0.125, 0.0, "limiter did not cap initial step");
  near(stepper.y_new[0], 0.125, 1e-15, "limited integration");
  require(stepper.stats.derivative_calls == 7 && stepper.stats.accepted_steps == 1 &&
              stepper.stats.rejected_steps == 0 && stepper.stats.limiter_calls == 1,
          "incorrect initial FSAL or limiter statistics");
  stepper.update_old();
  stepper.do_step();
  near(stepper.h_old, 0.125, 0.0, "limiter not reapplied after adaptive growth");
  require(stepper.stats.derivative_calls == 13, "FSAL derivative was not reused");
  stepper.init(1.0, -1.0, {1.0});
  stepper.do_step();
  near(stepper.h_old, -0.125, 0.0, "limiter changed integration direction");
  near(stepper.y_new[0], 0.875, 1e-15, "backward limited integration");
  stepper.reset_statistics();
  require(stepper.stats.derivative_calls == 0 && stepper.stats.accepted_steps == 0,
          "statistics reset failed");
  stepper.set_step_limiter({});
  stepper.init(0.0, 1.0, {0.0});
  stepper.do_step();
  near(stepper.h_old, 1.0, 0.0, "empty limiter should remove the bound");
}

void test_narrow_layer() {
  // Every stage of an unrestricted step from 0 to 1 misses this smooth layer.
  // Its exact integral is one: integral[-1,1] (1-u*u)^4 du = 256/315.
  constexpr double width = 0.01;
  auto layer = [](double x, State const &, State &dy) {
    double u = (x - 0.5) / width;
    double profile = std::max(0.0, 1.0 - u * u);
    dy[0] = 315.0 / (256.0 * width) * std::pow(profile, 4);
  };
  StepperDopr5<1, decltype(layer)> raw(layer, 1e-11, 1e-11);
  raw.init(0.0, 1.0, {0.0});
  raw.do_step();
  near(raw.y_new[0], 0.0, 0.0, "test must reproduce a completely skipped layer");
  for (double direction : {1.0, -1.0}) {
    StepperDopr5<1, decltype(layer)> limited(layer, 1e-11, 1e-11);
    const double end = direction > 0.0 ? 1.0 : 0.0;
    limited.set_step_limiter([=](double x, State const &) {
      return std::min(width / 4.0, std::abs(end - x));
    });
    limited.init(direction > 0.0 ? 0.0 : 1.0, direction, {0.0});
    unsigned iterations = 0;
    while (direction * (end - limited.x_old) > 0.0) {
      require(++iterations < 2000, "limited layer integration did not terminate");
      limited.do_step();
      limited.update_old();
    }
    near(limited.y_old[0], direction, 2e-8, "limited layer integral disagrees with analytic value");
    require(limited.stats.rejected_steps > 0, "layer should exercise rejected adaptive steps");
    require(limited.stats.derivative_calls ==
                1 + 6 * (limited.stats.accepted_steps + limited.stats.rejected_steps),
            "rejected steps missing from derivative statistics");
  }
}

void test_events_and_dense_output() {
  auto quadratic = [](double x, State const &, State &dy) { dy[0] = 2.0 * x; };
  StepperDopr5<1, decltype(quadratic)> stepper(quadratic, 1e-10, 1e-10);
  stepper.add_event([](double x, State const &) { return x - 0.8; });
  stepper.add_event([](double, State const &y) { return y[0] - 0.25; });
  stepper.add_event([](double x, State const &) { return x - 0.5; });
  stepper.init(0.0, 1.0, {0.0});
  stepper.do_step();
  require(stepper.detect_event() == 1, "earliest event was not selected");
  near(stepper.h_old, 0.5, 3e-15, "event location");
  near(stepper.y_new[0], 0.25, 3e-15, "truncated event state");
  near(stepper.dydx_new[0], 1.0, 6e-15, "event endpoint derivative");
  require(!stepper.events[0].active && stepper.events[1].active && stepper.events[2].active,
          "active flags must refer to the selected event location");
  stepper.prepare_dense();
  near(stepper.dense_out(0.25)[0], 0.0625, 3e-15,
       "dense output changed after event truncated the accepted step");
  stepper.update_old();
  near(stepper.dense_out(0.25)[0], 0.0625, 3e-15,
       "dense output lost its interval after update_old");
  stepper.h_old = 0.5;
  stepper.do_step();
  require(stepper.detect_event() == 0, "continuation retriggered the old event or lost a later event");
  near(stepper.x_old + stepper.h_old, 0.8, 3e-15, "continued event position");

  StepperDopr5<1, decltype(quadratic)> backward(quadratic, 1e-10, 1e-10);
  backward.init(1.0, -1.0, {1.0});
  // Adding events after initialization must initialize their signs too.
  backward.add_event([](double x, State const &) { return x - 0.4; });
  backward.add_event([](double x, State const &) { return x - 0.7; });
  backward.do_step();
  require(backward.detect_event() == 1, "backward event ordering is incorrect");
  near(backward.x_old + backward.h_old, 0.7, 3e-15, "backward event position");
  near(backward.dense_out(0.85)[0], 0.85 * 0.85, 3e-15, "backward truncated dense output");
}

void test_invalid_values() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  auto constant = [](double, State const &, State &dy) { dy[0] = 1.0; };
  throws([&] { StepperDopr5<1, decltype(constant)> bad(constant, 0.0, 0.0); },
         "zero tolerances accepted");
  StepperDopr5<1, decltype(constant)> stepper(constant, 1e-10, 1e-10);
  throws([&] { stepper.do_step(); }, "uninitialized step accepted");
  throws([&] { stepper.init(0.0, 0.0, {0.0}); }, "zero initial step accepted");
  throws([&] { stepper.init(nan, 1.0, {0.0}); }, "NaN initial abscissa accepted");
  throws([&] { stepper.init(0.0, 1.0, {nan}); }, "NaN initial state accepted");
  for (double bound : {0.0, -1.0, nan, -infinity}) {
    stepper.init(0.0, 1.0, {0.0});
    stepper.set_step_limiter([=](double, State const &) { return bound; });
    throws([&] { stepper.do_step(); }, "invalid limiter bound accepted");
  }
  stepper.set_step_limiter([=](double, State const &) { return infinity; });
  stepper.init(0.0, 1.0, {0.0});
  stepper.do_step();
  near(stepper.h_old, 1.0, 0.0, "positive infinite limiter should allow unrestricted steps");
  stepper.init(1e16, 0.01, {0.0});
  throws([&] { stepper.do_step(); }, "step that cannot advance floating point x accepted");
  auto invalid_derivative = [=](double x, State const &, State &dy) {
    dy[0] = x == 0.0 ? 0.0 : nan;
  };
  StepperDopr5<1, decltype(invalid_derivative)> bad_derivative(invalid_derivative, 1e-10, 1e-10);
  bad_derivative.init(0.0, 1.0, {0.0});
  throws([&] { bad_derivative.do_step(); }, "NaN stage derivative accepted");
  require(bad_derivative.stats.derivative_calls == 2,
          "non-finite stage derivative did not fail immediately");
  stepper.init(0.0, 1.0, {0.0});
  stepper.add_event([=](double x, State const &) { return x == 0.0 ? -1.0 : nan; });
  stepper.do_step();
  throws([&] { stepper.detect_event(); }, "NaN event accepted");
}

void test_roots() {
  auto scaled = [](double x) { return 1e300 * (x * x - 2.0); };
  near(zriddr(scaled, 2.0, 0.0, 1e-14), std::sqrt(2.0), 2e-14,
       "root with reversed interval and large function magnitude");
  near(zriddr([](double x) { return 1e-300 * (x - 0.3); }, 0.0, 1.0, 1e-14),
       0.3, 2e-14, "root with tiny function magnitude");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  throws([&] { zriddr([=](double) { return nan; }, 0.0, 1.0, 1e-8); },
         "NaN root function accepted");
  throws([&] { zriddr([](double x) { return x; }, -1.0, 1.0, -1.0); },
         "negative root tolerance accepted");
  throws([&] { zriddr([](double x) { return x * x + 1.0; }, -1.0, 1.0, 1e-8); },
         "unbracketed root accepted");
}
} // namespace

int main() {
  try {
    test_limiter_and_statistics();
    test_narrow_layer();
    test_events_and_dense_output();
    test_invalid_values();
    test_roots();
    std::cout << "Dopr5 regression tests passed\n";
    return 0;
  } catch (std::exception const &error) {
    std::cerr << "Dopr5 regression test failed: " << error.what() << '\n';
    return 1;
  }
}
