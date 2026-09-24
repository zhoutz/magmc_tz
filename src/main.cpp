#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "dopr5.hpp"
#include "photon.hpp"
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
};

double event_escape(double x, YVector const &y) {
  double r = y[0];
  return r - 10000;
}

BField bfield("table/bfield_t10.txt", B_pole, R_star);

int main() {
  for (double b0 : {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9}) {
    Boltzmann fb(-0.75);
    for (double muz : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}) {
      double3 r_hat{std::sqrt(1 - muz * muz), 0, muz};
      double3 n{0, 1, 0};
      for (double oi : {0.01, 0.1, 1., 10., 100.}) {
        for (Polarization pol : {Polarization::E, Polarization::O}) {
          PhotonEvolution photon_evolution{
              .bfield = bfield,
              .fb = fb,
              .n = n,
              .e1 = r_hat,
              .e2 = cross(n, r_hat),
              .r_psi_alpha_tau = {R_star, 0.0, 0.0, 0.0},
              .omega_inf = oi,
              .pol = pol,
          };

          StepperDopr5<4, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);
          stepper.add_event(&event_escape);
          photon_evolution.r_psi_alpha_tau[3] = 0;
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
                std::println("Photon escaped at r = {}, psi = {}, alpha = {}, tau = {}", r, psi,
                             alpha, tau);
                break;
              }
            }
            stepper.update_old();
          }
          double tau = stepper.y_new[3];

          std::println("b0={}, muz={}, pol={}, omega_inf={}, tau={}", b0, muz,
                       (pol == Polarization::E ? "E" : "O"), oi, tau);
        }
      }
    }
  }
}
