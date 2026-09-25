#pragma once

#include "bessel_K1_scaled.hpp"
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

struct Boltzmann {
  double b0, exp_k1, a;
  double b_min, b_max, b_mean;
  std::vector<double> knots;

  Boltzmann(double b0, int n_knots) : b0(b0) {
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

    knots.resize(n_knots);
    for (int i = 0; i < n_knots; ++i) {
      knots[i] = b_min + (b_max - b_min) * i / (n_knots - 1);
    }

    double s0 = std::sqrt((1 - b0) * (1 + b0));
    a = s0 * (1 + s0) / (b0 * b0);
    exp_k1 = bessel_K1_scaled(a);
    b_mean = std::copysign(1 / (a * exp_k1), b0);
  }

  double f(double b) const {
    if (!std::isfinite(b) || b <= b_min || b >= b_max) {
      return 0;
    }
    double s = std::sqrt((1 - b) * (1 + b));
    double gm1 = b * b / (s * (1 + s));
    return std::exp(-a * gm1) / (exp_k1 * s * s * s);
  }
};
