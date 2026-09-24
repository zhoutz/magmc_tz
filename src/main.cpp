#include "radial_optical_depth.hpp"

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
          // This executable launches outward radial rays in the northern
          // hemisphere. Seed the resonant layer and remove its endpoint
          // singularity before applying adaptive quadrature.
          RadialOpticalDepth integral(bfield, fb, muz, oi, pol);
          double tau = integral.integrate_to(10000, 1e-10, 1e-9);

          std::println("b0={}, muz={}, pol={}, omega_inf={}, tau={}", b0, muz,
                       (pol == Polarization::E ? "E" : "O"), oi, tau);
        }
      }
    }
  }
}
