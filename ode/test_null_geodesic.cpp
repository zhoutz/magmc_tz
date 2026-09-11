#include "dopr5.hpp"
#include <array>
#include <cmath>
#include <print>

constexpr int N = 3;
constexpr double rs = 1;
using YVector = std::array<double, N>;

void derivs(double x, YVector y, YVector &dydx) {
  auto [r, phi, alpha] = y;
  double f = std::sqrt(1 - rs / r);

  dydx[0] = f * std::cos(alpha);
  dydx[1] = std::sin(alpha) / r;
  dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
}

int main() {
  double x_init = 0.0;
  double h_init = 0.1;
  YVector y_init = {10.0, 0.0, M_PI / 4};
  YVector dydx_init;
  derivs(x_init, y_init, dydx_init);
  double atol = 1e-6;
  double rtol = 1e-6;

  StepperDopr5<N, decltype(derivs)> stepper(derivs, x_init, h_init, y_init, dydx_init, atol, rtol);

  for (int i = 0; i < 100; ++i) {
    stepper.do_step();
    stepper.update_old();
    auto [r, phi, alpha] = stepper.y_new;
    std::print("Step {}: r = {}, phi = {}, alpha = {}\n", i + 1, r, phi, alpha);
  }

  return 0;
}
