#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
#include "../qags.hpp"
#include "../solve_quadratic.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <print>
#include <random>

constexpr double M_star = 1.4;
constexpr double R_star = 10;
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km;
constexpr double B_pole = 1e14;

using YVector = std::array<double, 3>;

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  void operator()(double x, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha] = y;
    double f = std::sqrt(1 - rs / r);

    dydx[0] = f * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
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
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu_in =
        b.x * std::cos(alpha) + std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat));
    std::array<double, 2> betas;
    if (!solve_quadratic(x * x + mu_in * mu_in, -2 * mu_in, (1 + x) * (1 - x), betas)) return 0;
    double ret = 0;
    for (double beta : betas) {
      double f = fb.f(beta);
      if (f == 0) continue;
      double mup_in = (mu_in - beta) / (1 - beta * mu_in);
      double esq = (pol == Polarization::E) ? (0.5) : (0.5 * mup_in * mup_in);
      ret += f * esq * (1 - beta * mu_in) * (1 - beta * mu_in) * (1 + beta) * (1 - beta) /
             std::abs(mu_in - beta);
    }
    ret *= (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz) / (std::abs(fb.b_bar()) * r);

    if (!std::isfinite(ret)) {
      throw std::runtime_error("Non-finite dtaudl encountered");
    }

    return ret;
  }

  double calc_discriminant(double r, double psi, double alpha) const {
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = r_hat.z;
    double3 B_vec = bfield.calc_B(r, muz);
    double B = B_vec.length();
    double3 b = B_vec / B;
    double omega_c = B_to_omega * B;
    double omega = omega_inf / std::sqrt(1 - rs / r);
    double x = omega_c / omega;
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu_in =
        b.x * std::cos(alpha) + std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat));

    return x * x + mu_in * mu_in - 1;
  }
};

double event_escape(double x, YVector const &y) {
  double r = y[0];
  return r - 10000;
}

BField bfield("table/bfield_t10.txt", B_pole, R_star);

// Compute optical depth for non-radial trajectory
double total_optical_depth_nonradial(double b0, double initial_alpha, double muz, double oi, Polarization pol) {
  Boltzmann fb(b0);

  // Initial radial direction
  double3 r_hat{std::sqrt(1 - muz * muz), 0, muz};

  // Initial propagation direction (non-radial if alpha != 0)
  double cos_alpha = std::cos(initial_alpha);
  double sin_alpha = std::sin(initial_alpha);

  // Orbital plane normal: perpendicular to initial r and k
  // For simplicity, choose n in the x-z plane
  double3 n{-muz, 0, std::sqrt(1 - muz * muz)};
  n = n / n.length();

  double3 e1 = r_hat;
  double3 e2 = cross(n, r_hat);

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .n = n,
      .e1 = e1,
      .e2 = e2,
      .omega_inf = oi,
      .pol = pol,
  };

  StepperDopr5<3, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);
  stepper.add_event(&event_escape);

  stepper.add_event([&](double, YVector const &y) {
    auto [r, psi, alpha] = y;
    return photon_evolution.calc_discriminant(r, psi, alpha);
  });

  stepper.init(0.0, 1e-3 * R_star, {R_star, 0.0, initial_alpha});
  double tau = 0;

  while (true) {
    // Adaptive step limiting based on resonance proximity
    auto [r, psi, alpha] = stepper.y_old;
    double D = photon_evolution.calc_discriminant(r, psi, alpha);

    double max_step = 1e-2 * r;
    if (std::abs(D) < 0.1) {
      max_step = std::min(max_step, 1e-4 * r);
    } else if (std::abs(D) < 1.0) {
      max_step = std::min(max_step, 1e-3 * r);
    }

    stepper.do_step(max_step);
    stepper.prepare_dense();
    int event_id = stepper.detect_event();

    double dtau = qags(
        [&](double x) {
          auto [r, psi, alpha] = stepper.dense_out(x);
          return photon_evolution.calc_dtaudl(photon_evolution.n, photon_evolution.e1,
                                              photon_evolution.e2, r, psi, alpha,
                                              photon_evolution.omega_inf, photon_evolution.pol);
        },
        stepper.x_old, stepper.x_old + stepper.h_old, 1e-6, 1e-6);

    tau += dtau;

    if (event_id == 0 || event_id == 1) break;
    stepper.update_old();
  }
  return tau;
}

int main() {
  std::FILE *fp = std::fopen("output/test_nonradial_step_limiter.txt", "w");

  // Test non-radial trajectories with various initial angles
  std::vector<double> initial_alphas = {0.0, 0.1, 0.3, 0.5};  // radians

  std::println("Testing step_limiter method on non-radial photon paths...\n");

  for (double initial_alpha : initial_alphas) {
    for (double b0 : {-0.1, -0.3, -0.5}) {
      for (double muz : {0.0, 0.3, 0.6}) {
        for (double oi : {0.01, 1.0, 100.0}) {
          for (Polarization pol : {Polarization::E, Polarization::O}) {
            double tau = total_optical_depth_nonradial(b0, initial_alpha, muz, oi, pol);
            std::println(fp, "{:.2f} {:.2f} {:.2f} {:.2f} {} {:.16e}",
                        b0, initial_alpha, muz, oi,
                        (pol == Polarization::E ? 1 : 0), tau);
            std::println("alpha={:.2f}, b0={:.2f}, muz={:.2f}, oi={:.2f}, pol={}, tau={:.6e}",
                        initial_alpha, b0, muz, oi, (pol == Polarization::E ? "E" : "O"), tau);
          }
        }
      }
    }
  }

  std::fclose(fp);
  std::println("\nNon-radial test completed. Results saved to output/test_nonradial_step_limiter.txt");
}
