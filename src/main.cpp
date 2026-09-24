#include "transport.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <print>

BField bfield("table/bfield_t10.txt", B_pole, R_star);

int main() {
  for (double b0 : {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9}) {
    Boltzmann fb(b0);
    for (double muz : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}) {
      for (double oi : {0.01, 0.1, 1., 10., 100.}) {
        for (Polarization pol : {Polarization::E, Polarization::O}) {
          auto photon = transport::make_photon(R_star, muz, 0, 0, oi, pol);
          transport::Options options;
          options.method = transport::Method::event_guard;
          options.tolerance = 1e-6;
          // The same integrator accepts arbitrary alpha, azimuth, hemisphere,
          // initial radius and both outward and inward propagation.
          double tau = transport::integrate(bfield, fb, photon, options).tau;

          std::println("b0={}, muz={}, pol={}, omega_inf={}, tau={}", b0, muz,
                       (pol == Polarization::E ? "E" : "O"), oi, tau);
        }
      }
    }
  }
}
