#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "dopr5.hpp"
#include "init.hpp"
#include "photon.hpp"
#include "ran.hpp"
#include "solve_quadratic.hpp"

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
  Photon photon;
  Ran &ran;

  void operator()(double x, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha, tau] = y;
    double f = std::sqrt(1 - rs / r);

    dydx[0] = f * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
    dydx[3] = 0; // Assuming tau is constant or not evolving
  }

  double calc_dtaudl(double3 n, double3 e1, double3 e2, double r, double psi, double alpha,
                     double omega_inf, Polarization pol) {
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
    double beta1, beta2;
    if (!solve_quadratic(x * x + mu_in * mu_in, -2 * mu_in, 1 - x * x, beta1, beta2)) return 0;
    double ret = 0;
    for (double beta : {beta1, beta2}) {
      double f = fb.f(beta);
      double mu_in_p = (mu_in - beta) / (1 - beta * mu_in);
      double e2 = (pol == Polarization::E) ? (0.5) : (0.5 * mu_in_p * mu_in_p);
      ret += f * e2 * (1 - beta * mu_in) * (1 - beta * mu_in) * (1 - beta * beta) /
             std::abs(mu_in - beta);
    }
    ret *= (bfield.p + 1) * pi * b.z / (std::abs(fb.b_bar()) * r * b.y);
    return ret;
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

  StepperDopr5<4, PhotonEvolution> stepper(
      PhotonEvolution{.bfield = bfield, .fb = fb, .photon = photon, .ran = ran}, 1e-6, 1e-6);
  stepper.add_event(&event_escape);
  stepper.add_event(&event_absorption);
  stepper.add_event(&event_scattering);
  double target_tau = -std::log(ran.U());
  stepper.init(0.0, 1e-3 * R_star, YVector{photon.r, photon.psi, photon.alpha, -target_tau});

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
        break;
      }
    }
    stepper.update_old();
  }
}
