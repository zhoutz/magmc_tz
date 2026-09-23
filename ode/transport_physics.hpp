#pragma once

#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "photon.hpp"
#include "ran.hpp"
#include "sample_mup.hpp"
#include "solve_quadratic.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

constexpr double M_star = 1.4;                                    // M_sun
constexpr double R_star = 10;                                     // km
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km; // km
constexpr double B_pole = 1e14;                                   // G

using YVector = std::array<double, 4>;

struct ResonanceGeometry {
  double r, mu, x, discriminant, muz;
  double3 r_hat, b_cart;
};

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  Ran &ran;
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

  ResonanceGeometry geometry(YVector const &y) const {
    return geometry(n, e1, e2, y[0], y[1], y[2], omega_inf);
  }

  ResonanceGeometry geometry(double3 plane_n, double3 plane_e1, double3 plane_e2,
                            double r, double psi, double alpha, double energy) const {
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
    return {r, mu, x, std::fma(x, x, (mu - 1) * (mu + 1)), muz, r_hat, b_cart};
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

  void perform_scattering() {
    auto g = geometry(r_psi_alpha_tau);
    double3 r_hat = g.r_hat;
    double mu_in = g.mu;
    std::array<double, 2> betas{};
    auto weights = resonance_weights(g, pol, betas);
    double total_weight = weights[0] + weights[1];
    if (!(total_weight > 0.0) || !std::isfinite(total_weight)) {
      throw std::runtime_error("Invalid scattering weights");
    }
    double rand_val = ran.U() * total_weight;
    double beta = (rand_val < weights[0]) ? betas[0] : betas[1];

    double mup_out = sample_mup(ran);
    double mu_out = (mup_out + beta) / (1 + beta * mup_out);
    double3 b_cart = g.b_cart;
    double3 t_hat = ran.unit_perp_to(b_cart);
    double3 k_out = mu_out * b_cart + std::sqrt(1 - mu_out * mu_out) * t_hat;
    double alpha_out = std::atan2(cross(k_out, r_hat).length(), dot(k_out, r_hat));
    double3 e1_out = r_hat;
    double3 normal = cross(r_hat, k_out);
    double3 n_out = normal.length() > 1e-14 ? to_unit(normal) : ran.unit_perp_to(r_hat);
    double3 e2_out = cross(n_out, e1_out);
    double omega_inf_out = omega_inf * (1 - beta * mu_in) / (1 - beta * mu_out);
    Polarization pol_out =
        (ran.U() < 1 / (1 + mup_out * mup_out)) ? Polarization::E : Polarization::O;

    n = n_out;
    e1 = e1_out;
    e2 = e2_out;
    r_psi_alpha_tau[1] = 0;
    r_psi_alpha_tau[2] = alpha_out;
    r_psi_alpha_tau[3] = 0;
    pol = pol_out;
    omega_inf = omega_inf_out;

    if (!std::isfinite(omega_inf)) {
      throw std::runtime_error("Non-finite omega_inf encountered after scattering");
    }
  }
};
