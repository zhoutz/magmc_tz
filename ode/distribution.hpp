#pragma once

#include "bessel_K1_scaled.hpp"
#include <cmath>
#include <cstdlib>

struct Boltzmann {
  double b0, exp_k1, a, _b_bar;

  Boltzmann(double b0) : b0(b0) {
    double g0 = 1 / std::sqrt(1 - b0 * b0);
    a = 1 / (g0 - 1);
    exp_k1 = bessel_K1_scaled(a);
    _b_bar = std::copysign(1 / (a * exp_k1), b0);
  }

  double f(double b) const {
    if (b * b0 < 0) return 0;
    if (std::abs(b) >= 1) return 0;
    double g = 1 / std::sqrt(1 - b * b);
    double g3 = g * g * g;
    return std::exp(-a * (g - 1)) * g3 / exp_k1;
  }

  double b_bar() const { return _b_bar; }
};
