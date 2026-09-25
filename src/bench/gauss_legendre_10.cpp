#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
#include "../solve_quadratic.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <print>

constexpr double M_star = 1.4;
constexpr double R_star = 10;
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km;
constexpr double B_pole = 1e14;

using YVector = std::array<double, 3>;

// 10-point Gauss-Legendre quadrature nodes and weights on [-1, 1]
// Exact for polynomials up to degree 2*10-1 = 19
constexpr int GL_N = 10;
constexpr std::array<double, GL_N> GL_NODES = {
    -0.9739065285171717,
    -0.8650633666889845,
    -0.6794095682990244,
    -0.4333953941292472,
    -0.1488743389816312,
     0.1488743389816312,
     0.4333953941292472,
     0.6794095682990244,
     0.8650633666889845,
     0.9739065285171717
};
constexpr std::array<double, GL_N> GL_WEIGHTS = {
    0.0666713443086881,
    0.1494513491505806,
    0.2190863625159821,
    0.2692667193099964,
    0.2955242247147529,
    0.2955242247147529,
    0.2692667193099964,
    0.2190863625159821,
    0.1494513491505806,
    0.0666713443086881
};

template<typename F>
double gauss_legendre_10(F&& f, double a, double b) {
  double mid = 0.5 * (a + b);
  double half_len = 0.5 * (b - a);
  double sum = 0.0;
  for (int i = 0; i < GL_N; ++i) {
    double x = mid + half_len * GL_NODES[i];
    sum += GL_WEIGHTS[i] * f(x);
  }
  return half_len * sum;
}

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
};

double event_escape(double x, YVector const &y) {
  double r = y[0];
  return r - 10000;
}

BField bfield("table/bfield_t10.txt", B_pole, R_star);

double total_optical_depth(double b0, double muz, double oi, Polarization pol) {
  Boltzmann fb(b0);
  double3 r_hat{std::sqrt(1 - muz * muz), 0, muz};
  double3 n{0, 1, 0};

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .n = n,
      .e1 = r_hat,
      .e2 = cross(n, r_hat),
      .omega_inf = oi,
      .pol = pol,
  };

  StepperDopr5<3, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);
  stepper.add_event(&event_escape);

  stepper.add_event([&](double, YVector const &y) {
    auto [r, psi, alpha] = y;
    auto e1 = photon_evolution.e1;
    auto e2 = photon_evolution.e2;
    auto omega_inf = photon_evolution.omega_inf;

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

    return std::fma(x, x, (mu_in - 1) * (mu_in + 1));
  });

  stepper.add_event([&](double, YVector const &y) {
    double r = y[0];
    return B_to_omega * bfield.calc_B(r, muz).length() * std::sqrt(1 - rs / r) / oi - 1;
  });

  const bool stop_at_beta_zero = b0 * bfield.calc_B(R_star, muz).x <= 0;
  stepper.init(0.0, 1e-3 * R_star, {R_star, 0.0, 0.0});
  double tau = 0;

  while (true) {
    stepper.do_step(1e-2 * stepper.y_old[0]);
    stepper.prepare_dense();
    int event_id = stepper.detect_event();

    // Use 10-point Gauss-Legendre for higher accuracy
    double dtau = gauss_legendre_10(
        [&](double x) {
          auto [r, psi, alpha] = stepper.dense_out(x);
          return photon_evolution.calc_dtaudl(photon_evolution.n, photon_evolution.e1,
                                              photon_evolution.e2, r, psi, alpha,
                                              photon_evolution.omega_inf, photon_evolution.pol);
        },
        stepper.x_old, stepper.x_old + stepper.h_old);

    tau += dtau;

    if (event_id == 0 || event_id == 1 || (event_id == 2 && stop_at_beta_zero)) break;
    stepper.update_old();
  }
  return tau;
};

int main() {
  std::FILE *fp = std::fopen("output/gauss_legendre_10.txt", "w");
  for (double b0 : {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9}) {
    for (double muz : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}) {
      for (double oi : {0.01, 0.1, 1., 10., 100.}) {
        for (Polarization pol : {Polarization::E, Polarization::O}) {
          double tau = total_optical_depth(b0, muz, oi, pol);
          std::println(fp, "{:.2f} {:.2f} {:.2f} {} {:.16e}", b0, muz, oi,
                       (pol == Polarization::E ? 1 : 0), tau);
        }
      }
    }
  }
  std::fclose(fp);
}
