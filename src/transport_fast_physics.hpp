#pragma once
// Deterministic physics subset ported from codex-run1 commit
// f373eb507956aeb22eff5e16577573d26b171efc, ode/transport_physics.hpp.
// Only the unused RNG member and perform_scattering() were removed.
#include "photon_evolution.hpp"
#include <algorithm>
#include <stdexcept>
namespace fast_transport {
struct ResonanceGeometry {
  double r, mu, x, discriminant, muz;
  double3 r_hat, b_cart;
  double current = 0; // Optional cached Bphi/Btheta, used by fast transport.
};

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  double3 n, e1, e2;
  YVector r_psi_alpha_tau;
  double omega_inf;
  Polarization pol;

  void operator()(double, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha, tau] = y;
    double f = std::sqrt(1 - rs / r);

    dydx[0] = f * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
    dydx[3] = calc_dtaudl(n, e1, e2, r, psi, alpha, omega_inf, pol);
  }

  ResonanceGeometry geometry(YVector const &y, bool cache_current = false) const {
    return geometry(n, e1, e2, y[0], y[1], y[2], omega_inf, cache_current);
  }

  ResonanceGeometry geometry(double3 plane_n, double3 plane_e1, double3 plane_e2,
                            double r, double psi, double alpha, double energy,
                            bool cache_current = false) const {
    if (!(r > rs) || !(energy > 0) || !std::isfinite(energy))
      throw std::runtime_error("Invalid photon radius or energy");
    double3 r_hat = to_unit(std::cos(psi) * plane_e1 + std::sin(psi) * plane_e2);
    double muz = std::clamp(r_hat.z, -1.0, 1.0);
    double3 B_sph = bfield.calc_B(r, muz);
    double B = B_sph.length();
    double rho = std::hypot(r_hat.x, r_hat.y);
    // The transverse field vanishes on the axis, so any azimuthal basis works.
    double3 theta_hat = rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                               : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double3 b_cart = (B_sph.x * r_hat + B_sph.y * theta_hat + B_sph.z * phi_hat) / B;
    double3 k = std::cos(alpha) * r_hat + std::sin(alpha) * cross(plane_n, r_hat);
    double mu = std::clamp(dot(k, b_cart), -1.0, 1.0);
    double x = B_to_omega * B * std::sqrt(1 - rs / r) / energy;
    double current = cache_current && B_sph.y != 0 ? B_sph.z / B_sph.y : 0;
    return {r, mu, x, std::fma(x, x, (mu - 1) * (mu + 1)), muz, r_hat, b_cart, current};
  }

  std::array<double, 2> resonance_weights(ResonanceGeometry const &g,
                                        Polarization mode,
                                        std::array<double, 2> &betas) const {
    std::array<double, 2> weights{};
    if (!(g.discriminant > 0)) return weights; // Exact coalescence is integrated by a change of variable.
    double sd = std::sqrt(g.discriminant);
    double aa = g.x * g.x + g.mu * g.mu;
    double q = g.mu + std::copysign(g.x * sd, g.mu);
    if (q == 0) return weights;
    betas = {q / aa, std::abs(g.mu) < .5 * g.x * sd
                        ? (g.mu - std::copysign(g.x * sd, g.mu)) / aa
                        : ((1 - g.x) * (1 + g.x)) / q};
    for (int i = 0; i < 2; ++i) {
      double beta = betas[i];
      if (!std::isfinite(beta) || !(std::abs(beta) < 1)) continue;
      double f = fb.f(beta);
      if (f == 0) continue;
      // |mu'| = sqrt(D)/x avoids subtracting nearly coalescing roots.
      weights[i] = 0.5 * f * (1 - beta * g.mu) * (1 - beta * beta) *
                   (mode == Polarization::E ? g.x / sd : sd / g.x);
    }
    return weights;
  }

  double rate(ResonanceGeometry const &g, Polarization mode) const {
    double current = bfield.Bphi_over_Btheta(g.muz);
    if (current == 0) return 0;
    std::array<double, 2> betas{};
    auto weights = resonance_weights(g, mode, betas);
    double ret = (weights[0] + weights[1]) * (bfield.p + 1) * pi * current /
                 (std::abs(fb.b_bar()) * g.r);
    if (!(ret >= 0) || !std::isfinite(ret)) throw std::runtime_error("Non-finite dtaudl encountered");
    return ret;
  }

  double rate(ResonanceGeometry const &g) const { return rate(g, pol); }

  double calc_dtaudl(double3 plane_n, double3 plane_e1, double3 plane_e2, double r,
                     double psi, double alpha, double energy, Polarization mode) const {
    return rate(geometry(plane_n, plane_e1, plane_e2, r, psi, alpha, energy), mode);
  }

};
}
