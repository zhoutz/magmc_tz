#include "resonance_reference.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
  BField field("table/bfield_t10.txt", 1e14, 10);
  const double rs = 1.4 * schwarzschild_radius_of_sun_in_km;
  std::cout << std::setprecision(16);
  for (double b0 : {-0.75, -0.2, -0.05, -0.005, -0.001, 0.2}) {
    Boltzmann distribution(b0);
    for (double muz : {-0.2, 0.0}) {
      for (Polarization pol : {Polarization::O, Polarization::E}) {
        Photon photon{.n = {0, 1, 0},
                      .e1 = {std::sqrt(1 - muz * muz), 0, muz},
                      .e2 = {muz, 0, -std::sqrt(1 - muz * muz)},
                      .r = 10,
                      .psi = 0,
                      .alpha = 0,
                      .omega_inf = 1,
                      .pol = pol};
        const auto reference =
            resonance_reference::radial(field, distribution, rs, photon, 10000, 1e-10);
        const auto tighter =
            resonance_reference::radial(field, distribution, rs, photon, 10000, 1e-12);
        // Reverse direction and carrier velocities together: the radial resonance
        // and all overlaps are unchanged, hence so is total optical depth.
        Photon reversed = photon;
        reversed.r = 10000;
        reversed.alpha = pi;
        const auto reverse =
            resonance_reference::radial(field, Boltzmann(-b0), rs, reversed, 10, 1e-11);
        const double scale = std::max(1.0, std::abs(reference.tau));
        if (std::abs(reference.tau - tighter.tau) > 1e-10 * scale ||
            std::abs(reverse.tau - reference.tau) > 1e-10 * scale)
          throw std::runtime_error("independent radial reference consistency failure");
        // Additivity exercises support intervals clipped by propagation endpoints.
        photon.r = 10;
        const auto inner = resonance_reference::radial(field, distribution, rs, photon, 100, 1e-11);
        photon.r = 100;
        const auto outer =
            resonance_reference::radial(field, distribution, rs, photon, 10000, 1e-11);
        if (std::abs(inner.tau + outer.tau - reference.tau) > 1e-10 * scale)
          throw std::runtime_error("radial reference additivity failure");
        if (b0 == -0.2 && muz == -0.2) {
          const double expected = pol == Polarization::O ? 0.461715239204475 : 4.66235008107873;
          if (std::abs(tighter.tau - expected) > 2e-12)
            throw std::runtime_error("radial reference disagrees with independent review value");
        }
        if (muz == 0 && pol == Polarization::E) {
          photon.r = 10;
          const auto flat =
              resonance_reference::radial(field, distribution, 0, photon, 10000, 1e-11);
          const double exact = (field.p + 1) * pi * field.Bphi_over_Btheta(0) /
                               (2 * (2 + field.p) * std::abs(distribution.b_bar()));
          if (std::abs(flat.tau - exact) > 2e-12 * std::max(1.0, exact))
            throw std::runtime_error(
                "radial reference disagrees with flat equatorial E analytic value");
        }
        std::cout << b0 << " muz=" << muz << ' ' << (pol == Polarization::O ? 'O' : 'E') << ' '
                  << tighter.tau << " estimate=" << tighter.error
                  << " evaluations=" << tighter.evaluations << '\n';
      }
    }
  }
}
