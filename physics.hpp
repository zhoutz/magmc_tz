#pragma once

#include <algorithm>
#include <cmath>

#include "solve_quadratic.hpp"

inline double muval(double muk, double phk, double muz, double phi, double bmod, double br,
                    double bth, double bph) {
  const double sinthk = std::sqrt(std::abs(1.0 - muk * muk));
  const double sinthz = std::sqrt(std::abs(1.0 - muz * muz));
  const double aux1 = sinthk * sinthz * std::cos(phk - phi) + muk * muz;
  const double aux2 = sinthk * muz * std::cos(phk - phi) - muk * sinthz;
  const double aux3 = sinthk * std::sin(phk - phi);
  double mu = (br / bmod) * aux1 + (bth / bmod) * aux2 + (bph / bmod) * aux3;
  return mu;
}

inline auto gambrs(double bmod, double omega, double mu) {
  double a = bmod / omega;
  double A = 1 - mu * mu;
  double B = -2 * a * mu;
  double C = 1 - a * a;
  return solve_quadratic_real(A, B, C);
}

inline double disfpw(double gbres, double gbmin, double gbmax, bool unidir) {
  if (!unidir) gbres = std::abs(gbres);
  if (gbmin < gbres && gbres < gbmax) {
    double norm = gbmax * gbmin / (gbmax - gbmin);
    return (unidir ? 1.0 : 0.5) * norm / (gbres * gbres);
  }
  return 0;
}
