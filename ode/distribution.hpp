#pragma once

#include <cmath>
#include <cstdlib>

struct Boltzmann {
  double b0, g0, k1, a, _b_bar;

  Boltzmann(double b0) : b0(b0) {
    g0 = 1 / std::sqrt(1 - b0 * b0);
    a = 1 / (g0 - 1);
    k1 = std::cyl_bessel_k(1, a);
    _b_bar = std::exp(-a) / (a * k1);
    _b_bar = std::copysign(_b_bar, b0);
  }

  double f(double b) const {
    if (b * b0 < 0) return 0;
    if (std::abs(b) >= 1) return 0;
    double g = 1 / std::sqrt(1 - b * b);
    double g3 = g * g * g;
    return std::exp(-a * g) * g3 / k1;
  }

  double b_bar() const { return _b_bar; }
};
