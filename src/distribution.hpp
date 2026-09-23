#pragma once

#include "bessel_K1_scaled.hpp"
#include <cmath>
#include <cstdlib>
#include <stdexcept>

struct Boltzmann {
  double b0, exp_k1, a, _b_bar;
  double b_min, b_max;

  Boltzmann(double b0) : b0(b0) {
    if (!std::isfinite(b0) || std::abs(b0) >= 1 || b0 == 0) {
      throw std::runtime_error("Invalid b0 value for Boltzmann distribution");
    }
    if (b0 < 0) {
      b_min = -1;
      b_max = 0;
    } else {
      b_min = 0;
      b_max = 1;
    }

    double s0 = std::sqrt((1 - b0) * (1 + b0));
    a = s0 * (1 + s0) / (b0 * b0);
    exp_k1 = bessel_K1_scaled(a);
    _b_bar = std::copysign(1 / (a * exp_k1), b0);
  }

  double f(double b) const {
    if (b * b0 < 0) return 0;
    if (std::abs(b) >= 1) return 0;
    double s = std::sqrt((1 - b) * (1 + b));
    double gm1 = b * b / (s * (1 + s));
    return std::exp(-a * gm1) / (exp_k1 * s * s * s);
  }

  double b_bar() const { return _b_bar; }
};
