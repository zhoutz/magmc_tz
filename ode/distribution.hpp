#pragma once

#include <cmath>
#include <cstdlib>

struct Boltzmann {
  double b0, g0, k1;

  Boltzmann(double b0) : b0(b0) {
    g0 = 1 / std::sqrt(1 - b0 * b0);
    k1 = std::cyl_bessel_k(1, 1 / (g0 - 1));
  }

  double f(double b) const {
    if (b * b0 < 0) return 0;
    if (std::abs(b) >= 1) return 0;
    double t1 = std::sqrt(1 - b * b);
    double t2 = t1 * (1 - b * b);
    double g = 1 / t1;
    return std::exp(-g / (g0 - 1)) / (k1 * t2);
  }

  double b_bar() const { return b0; }
};
