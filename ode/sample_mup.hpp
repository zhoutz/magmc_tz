#pragma once

#include "ran.hpp"

#include <algorithm>
#include <cmath>

inline double sample_mup(Ran &ran) {
  const double u = ran.U();
  if (u < 0.75) return (8.0 / 3.0) * u - 1.0;
  const double y = 8.0 * u - 7.0;
  double m = std::max(std::abs(y), ran.U());
  m = std::max(m, ran.U());
  return std::copysign(m, y);
}
