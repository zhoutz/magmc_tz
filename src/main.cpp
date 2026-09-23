#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "dopr5.hpp"
#include "init.hpp"
#include "photon.hpp"
#include "ran.hpp"
#include "sample_mup.hpp"
#include "solve_quadratic.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <print>

constexpr double M_star = 1.4;                                    // M_sun
constexpr double R_star = 10;                                     // km
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km; // km
constexpr double B_pole = 1e14;                                   // G

using YVector = std::array<double, 4>;

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  Ran &ran;
  double3 n, e1, e2;
  YVector r_psi_alpha_tau;
  double omega_inf;
  Polarization pol;

  void operator()(double x, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha, tau] = y;
    double f = std::sqrt(1 - rs / r);

    dydx[0] = f * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
    dydx[3] = calc_dtaudl(n, e1, e2, r, psi, alpha, omega_inf, pol);
  }

  double calc_dtaudl(double3 n, double3 e1, double3 e2, double r, double psi, double alpha,
                     double omega_inf, Polarization pol) const {
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = r_hat.z;
    double3 B_vec = bfield.calc_B(r, muz);
    double B = B_vec.length();
    double3 b = B_vec / B;
    double omega_c = B_to_omega * B;
    double omega = omega_inf / std::sqrt(1 - rs / r);
    double x = omega_c / omega;
    double rho = std::sqrt(r_hat.x * r_hat.x + r_hat.y * r_hat.y);
    double3 theta_hat{r_hat.x * r_hat.z / rho, r_hat.y * r_hat.z / rho, -rho};
    double3 phi_hat{-r_hat.y / rho, r_hat.x / rho, 0};
    double mu_in =
        b.x * std::cos(alpha) + std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat));
    std::array<double, 2> betas;
    if (!solve_quadratic(x * x + mu_in * mu_in, -2 * mu_in, 1 - x * x, betas)) return 0;
    double ret = 0;
    for (double beta : betas) {
      double f = fb.f(beta);
      if (f == 0) continue;
      double mup_in = (mu_in - beta) / (1 - beta * mu_in);
      double esq = (pol == Polarization::E) ? (0.5) : (0.5 * mup_in * mup_in);
      ret += f * esq * (1 - beta * mu_in) * (1 - beta * mu_in) * (1 - beta * beta) /
             std::abs(mu_in - beta);
    }
    ret *= (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz) / (std::abs(fb.b_bar()) * r);

    if (!std::isfinite(ret)) {
      throw std::runtime_error("Non-finite dtaudl encountered");
    }

    return ret;
  }

  void perform_scattering() {
    double r = r_psi_alpha_tau[0];
    double psi = r_psi_alpha_tau[1];
    double alpha = r_psi_alpha_tau[2];

    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = r_hat.z;
    double3 B_sph = bfield.calc_B(r, muz);
    double B = B_sph.length();
    double3 b_sph = B_sph / B;
    double omega_c = B_to_omega * B;
    double omega = omega_inf / std::sqrt(1 - rs / r);
    double x = omega_c / omega;
    double rho = std::sqrt(r_hat.x * r_hat.x + r_hat.y * r_hat.y);
    double3 theta_hat{r_hat.x * r_hat.z / rho, r_hat.y * r_hat.z / rho, -rho};
    double3 phi_hat{-r_hat.y / rho, r_hat.x / rho, 0};
    double mu_in = b_sph.x * std::cos(alpha) +
                   std::sin(alpha) * (b_sph.y * dot(n, phi_hat) - b_sph.z * dot(n, theta_hat));
    std::array<double, 2> betas;
    if (!solve_quadratic(x * x + mu_in * mu_in, -2 * mu_in, 1 - x * x, betas)) {
      throw std::runtime_error("No valid beta found for scattering");
    }
    std::array<double, 2> weights{};
    for (int i = 0; i < 2; ++i) {
      double beta = betas[i];
      double f = fb.f(beta);
      if (f == 0) continue;
      double mu_in_p = (mu_in - beta) / (1 - beta * mu_in);
      double esq = (pol == Polarization::E) ? (0.5) : (0.5 * mu_in_p * mu_in_p);
      weights[i] = f * esq * (1 - beta * mu_in) * (1 - beta * mu_in) * (1 - beta * beta) /
                   std::abs(mu_in - beta);
    }
    double total_weight = weights[0] + weights[1];
    if (!(total_weight > 0.0) || !std::isfinite(total_weight)) {
      throw std::runtime_error("Invalid scattering weights");
    }
    double rand_val = ran.U() * total_weight;
    double beta = (rand_val < weights[0]) ? betas[0] : betas[1];

    double mup_out = sample_mup(ran);
    double mu_out = (mup_out + beta) / (1 + beta * mup_out);
    double3 b_cart = b_sph.x * r_hat + b_sph.y * theta_hat + b_sph.z * phi_hat;
    double3 t_hat = ran.unit_perp_to(b_cart);
    double3 k_out = mu_out * b_cart + std::sqrt(1 - mu_out * mu_out) * t_hat;
    double alpha_out = std::atan2(cross(k_out, r_hat).length(), dot(k_out, r_hat));
    double3 e1_out = r_hat;
    double3 n_out = to_unit(cross(r_hat, k_out));
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

double event_escape(double x, YVector const &y) {
  double r = y[0];
  return r - 1000 * R_star;
}

double event_absorption(double x, YVector const &y) {
  double r = y[0];
  return r - R_star;
}

double event_scattering(double x, YVector const &y) {
  double tau = y[3];
  return tau;
}

BField bfield("table/bfield_t10.txt", B_pole, R_star);
Boltzmann fb(-0.75);

int main() {
  Ran ran(1234);
  Photon photon = init07(ran, R_star, Polarization::O);
  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .ran = ran,
      .n = photon.n,
      .e1 = photon.e1,
      .e2 = photon.e2,
      .r_psi_alpha_tau = {photon.r, photon.psi, photon.alpha, 0.0},
      .omega_inf = photon.omega_inf,
      .pol = photon.pol,
  };

  StepperDopr5<4, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);
  stepper.add_event(&event_escape);
  stepper.add_event(&event_absorption);
  stepper.add_event(&event_scattering);
  double target_tau = std::log(ran.U());
  photon_evolution.r_psi_alpha_tau[3] = target_tau;
  stepper.init(0.0, 1e-3 * R_star, photon_evolution.r_psi_alpha_tau);

  while (true) {
    stepper.do_step();
    int event_id = stepper.detect_event();
    if (event_id != -1) {
      double r = stepper.y_new[0];
      double psi = stepper.y_new[1];
      double alpha = stepper.y_new[2];
      double tau = stepper.y_new[3];

      if (event_id == 0) {
        std::println("Photon escaped at r = {}, psi = {}, alpha = {}, tau = {}", r, psi, alpha,
                     tau);
        break;
      } else if (event_id == 1) {
        std::println("Photon absorbed at r = {}, psi = {}, alpha = {}, tau = {}", r, psi, alpha,
                     tau);
        break;
      } else if (event_id == 2) {
        std::println("Photon scattered at r = {}, psi = {}, alpha = {}, tau = {}", r, psi, alpha,
                     tau);
        photon_evolution.r_psi_alpha_tau = stepper.y_new;
        photon_evolution.perform_scattering();
        photon_evolution.r_psi_alpha_tau[3] = std::log(ran.U());
        stepper.init(0, stepper.h_new, photon_evolution.r_psi_alpha_tau);
        continue;
      }
    }
    stepper.update_old();
  }
}
