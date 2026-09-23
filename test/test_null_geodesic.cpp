#include "resonance_transport.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using State = resonance_transport_detail::GeometryState;

void require(bool condition, std::string const &message) {
  if (!condition) throw std::runtime_error(message);
}

double impact_parameter(State const &state) {
  return state[0] * std::sin(state[2]) / std::sqrt(1 - rs / state[0]);
}

double radial_distance_primitive(double r) {
  return std::sqrt(r * (r - rs)) + rs * std::log(std::sqrt(r) + std::sqrt(r - rs));
}

void check_ray(std::string const &name, State const &initial, bool absorbed, bool turns) {
  constexpr double outer_radius = 100 * R_star;
  TransportStats statistics;
  // Exercise the same geometry RHS used by production transport, without any
  // optical-depth kernel, field table, carrier distribution, or random draws.
  resonance_transport_detail::GeometryRHS geometry{statistics};
  StepperDopr5<3, decltype(geometry)> stepper(geometry, 1e-12, 1e-11);
  stepper.set_step_limiter([](double, State const &y) { return .05 * y[0]; });
  stepper.add_event([](double, State const &y) { return y[0] - outer_radius; });
  stepper.add_event([](double, State const &y) { return y[0] - R_star; });
  stepper.add_event([](double, State const &y) { return y[2] - pi / 2; });
  stepper.init(0, .01 * R_star, initial);
  double const initial_impact = impact_parameter(initial);
  double const impact_scale = std::max(R_star, std::abs(initial_impact));
  double maximum_error = 0;
  bool reached_boundary = false;
  unsigned turning_events = 0;
  auto check_state = [&](State const &y) {
    for (double value : y) require(std::isfinite(value), name + ": non-finite state");
    maximum_error = std::max(maximum_error,
                             std::abs(impact_parameter(y) - initial_impact) / impact_scale);
  };
  for (int iteration = 0; iteration < 10000; ++iteration) {
    stepper.do_step();
    int const event = stepper.detect_event();
    stepper.prepare_dense();
    check_state(stepper.y_new);
    check_state(stepper.dense_out(stepper.x_old + .5 * stepper.h_old));
    if (event == 2) {
      ++turning_events;
      require(turns, name + ": unexpected radial turning point");
      require(std::abs(stepper.y_new[2] - pi / 2) < 1e-10,
              name + ": turning event angle is inaccurate");
      // Independently solve r^3-b^2*r+b^2*rs=0 for the outer turning radius.
      double lo = std::max(R_star, 1.5 * rs), hi = initial[0];
      double const b2 = initial_impact * initial_impact;
      for (int i = 0; i < 80; ++i) {
        double const middle = .5 * (lo + hi);
        double const polynomial = middle * middle * middle - b2 * middle + b2 * rs;
        if (polynomial > 0) hi = middle; else lo = middle;
      }
      require(std::abs(stepper.y_new[0] - .5 * (lo + hi)) < 2e-7,
              name + ": turning radius disagrees with conserved impact parameter");
    } else if (event >= 0) {
      require(event == (absorbed ? 1 : 0), name + ": wrong boundary event");
      double const expected_radius = absorbed ? R_star : outer_radius;
      require(std::abs(stepper.y_new[0] - expected_radius) < 1e-8,
              name + ": boundary radius is inaccurate");
      require((turning_events == 1) == turns, name + ": wrong number of turning events");
      if (std::abs(std::sin(initial[2])) < 1e-14) {
        double const exact_distance = std::abs(radial_distance_primitive(expected_radius) -
                                                radial_distance_primitive(initial[0]));
        require(std::abs(stepper.x_old + stepper.h_old - exact_distance) < 2e-7,
                name + ": radial affine distance disagrees with analytic solution");
        require(std::abs(stepper.y_new[1] - initial[1]) < 1e-12,
                name + ": radial ray changed azimuth");
      }
      reached_boundary = true;
      break;
    }
    stepper.update_old();
  }
  require(reached_boundary, name + ": failed to reach boundary");
  require(maximum_error < 2e-8, name + ": Schwarzschild impact parameter drift is too large");
}
} // namespace

int main() {
  try {
    check_ray("radial outward", {R_star, .3, 0}, false, false);
    check_ray("radial inward", {3 * R_star, -.2, pi}, true, false);
    check_ray("oblique outward", {R_star, .1, 1.1}, false, false);
    check_ray("inward then turning", {6 * R_star, 0, 2.3}, false, true);
    check_ray("oblique absorption", {2 * R_star, .1, 3.0}, true, false);
    std::cout << "Null geodesic regression tests passed\n";
    return 0;
  } catch (std::exception const &error) {
    std::cerr << "Null geodesic regression test failed: " << error.what() << '\n';
    return 1;
  }
}
