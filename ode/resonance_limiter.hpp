#pragma once

#include "transport_physics.hpp"
#include <algorithm>
#include <cmath>

// A low-cost guard for ordinary coupled (r,psi,alpha,tau) RK integration.
// It resolves the *frequency* width of a thermal population even when both
// resonant roots currently have negligible density. It is not a certified
// enclosure of a nonmonotonic ray; event quadrature handles the E-mode fold.
struct ThermalStepLimiter {
  PhotonEvolution const &physics;
  double resolution = 0.1;
  double max_fraction = 0.1;

  double operator()(double, YVector const &y) const {
    auto g = physics.geometry(y);
    double eps = 1e-5 * y[0];
    double lapse = std::sqrt(1 - rs / y[0]);
    YVector probe = y;
    probe[0] += eps * lapse * std::cos(y[2]);
    probe[1] += eps * std::sin(y[2]) / y[0];
    probe[2] -= eps * std::sin(y[2]) * (1 - 1.5 * rs / y[0]) / (y[0] * lapse);
    auto gp = physics.geometry(probe);
    double dlogx = (std::log(gp.x) - std::log(g.x)) / eps;
    double dmu = (gp.mu - g.mu) / eps;
    auto beta_at = [&](double u) {
      double gm1 = u * u / physics.fb.a;
      return std::copysign(std::sqrt(gm1 * (2 + gm1)) / (1 + gm1), physics.fb.b0);
    };
    auto frequency = [&](double b) {
      return std::log1p(-b * g.mu) - 0.5 * std::log1p(-b * b);
    };
    double b1 = beta_at(1), b6 = beta_at(6);
    double f1 = frequency(b1), fhalf = frequency(beta_at(0.5));
    double width = std::max({0.0, f1, fhalf}) - std::min({0.0, f1, fhalf});
    double lower = std::min(0.0, frequency(b6));
    double upper = std::max(0.0, frequency(b6));
    if (g.mu * b6 > 0 && std::abs(g.mu) < std::abs(b6))
      lower = std::min(lower, frequency(g.mu));
    double logx = std::log(g.x);
    double distance = std::max({0.0, lower - logx, logx - upper});
    double speed = std::abs(dlogx) + std::abs(b6 * dmu) / std::max(1e-12, 1 - b6 * g.mu);
    double limit = resolution * (distance + width) / std::max(speed, 1e-15 / y[0]);
    return std::min(max_fraction * y[0], limit);
  }
};
