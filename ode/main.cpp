#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "dopr5.hpp"
#include "init.hpp"
#include "ran.hpp"

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
