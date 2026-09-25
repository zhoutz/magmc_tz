#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../quad.hpp"
#include "../ran.hpp"
#include "../sample_mup.hpp"
#include "../solve_quadratic.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <print>

constexpr double M_star = 1.4;                                    // M_sun
constexpr double R_star = 10;                                     // km
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km; // km
constexpr double B_pole = 1e14;                                   // G

using YVector = std::array<double, 3>;

enum class Polarization { O, E };

struct PhotonEvolution {
  static void geodesic(double, YVector const &y, YVector &dydx) {
    auto [r, psi, alpha] = y;
    double L = std::sqrt(1 - rs / r);
    dydx[0] = L * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * L) * (1 - 3 * rs / (2 * r));
  }
  static double event_absorb(double, YVector const &y) { return y[0] - R_star; }
  static double event_escape(double, YVector const &y) { return y[0] - 10000; }
  static auto compute_knots(Boltzmann const &fb, int n_knots) {
    std::vector<double> knots(n_knots);
    for (int i = 0; i < n_knots; ++i) {
      knots[i] = fb.b_min + (fb.b_max - fb.b_min) * i / (n_knots - 1);
    }
    return knots;
  }

  BField const &bfield;
  Boltzmann const &fb;
  StepperDopr5<3, decltype(geodesic)> stepper;
  Ran ran;

  std::vector<double> knots;
  const int n_scan;
  const double orbit_atol, orbit_rtol, quad_atol, quad_rtol;

  Quad quad;

  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  PhotonEvolution(BField const &bfield, Boltzmann const &fb, int seed, //
                  int n_knots = 4, int n_scan = 4,                     //
                  double orbit_atol = 1e-10, double orbit_rtol = 1e-10,
                  double quad_atol = 1e-8, double quad_rtol = 1e-8)
      : bfield(bfield), fb(fb), stepper(geodesic, orbit_atol, orbit_rtol),
        ran(seed), knots(compute_knots(fb, n_knots)), n_scan(n_scan),
        orbit_atol(orbit_atol), orbit_rtol(orbit_rtol), quad_atol(quad_atol),
        quad_rtol(quad_rtol) {
    stepper.add_event(&event_absorb);
    stepper.add_event(&event_escape);
  }

  void init(double3 _n, double3 _e1, double3 _e2, double _omega_inf,
            Polarization _pol, double r0, double psi0, double alpha0) {
    n = _n;
    e1 = _e1;
    e2 = _e2;
    omega_inf = _omega_inf;
    pol = _pol;
    stepper.init(0.0, 1e-3 * R_star, {r0, psi0, alpha0});
  }

  void init_random_radial(double omega_inf, Polarization pol) {
    double3 e1 = ran.point_on_unit_sphere();
    double3 n = ran.unit_perp_to(e1);
    double3 e2 = cross(n, e1);
    init(n, e1, e2, omega_inf, pol, R_star, 0., 0.);
  }

  struct Geometry {
    double x, mu, D, pref;
  };

  Geometry geo(YVector const &y) const {
    auto [r, psi, alpha] = y;
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = std::clamp(r_hat.z, -1.0, 1.0);
    double3 B_vec = bfield.calc_B(r, muz);
    double B = B_vec.length();
    double3 b = B_vec / B;
    double x = B_to_omega * B * std::sqrt(1 - rs / r) / omega_inf;
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * r_hat.z / rho,
                                          r_hat.y * r_hat.z / rho, -rho}
                                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu = std::clamp(
        b.x * std::cos(alpha) +
            std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat)),
        -1.0, 1.0);
    double D = std::fma(x, x, (mu - 1) * (mu + 1));
    double overlap = (pol == Polarization::E) ? (0.5) : (0.5 * D / (x * x));
    double pref = (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz) *
                  fb.inv_abs_b_mean * x * x * overlap / (r * std::sqrt(D));
    if (!std::isfinite(pref) || pref < 0) {
      pref = 0;
    }
    return Geometry{
        .x = x,
        .mu = mu,
        .D = D,
        .pref = pref,
    };
  }

  auto rate(std::array<double, 2> const &betas) const {
    std::array<double, 2> ret{};
    for (int i = 0; i < 2; ++i) {
      double beta = betas[i];
      double f = fb.f(beta);
      if (f == 0)
        continue;
      double t1 = (1 + beta) * (1 - beta);
      ret[i] = f * t1 * std::sqrt(t1);
    }
    return ret;
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
  } escaped_data;

  struct ScatteredData {
    YVector r_psi_alpha;
    double beta;
  } scattered_data;

  EvolveResult evolve_geodesic() {
    double tau = std::log(ran.U());
    std::vector<double> cuts;

    while (true) {
      stepper.do_step(0.1 * stepper.y_old[0]);
      int event_id = stepper.detect_event();

      double xl = stepper.x_old, xr = stepper.x_old + stepper.h_old;
      cuts.clear(), cuts.push_back(xl), cuts.push_back(xr);
      double l = xl;
      auto gl = geo(stepper.dense_out(l));
      for (int i = 1; i <= n_scan; ++i) {
        double r = xl + (xr - xl) * i / n_scan;
        auto gr = geo(stepper.dense_out(r));
        if (gl.D * gr.D <= 0) {
          cuts.push_back(zriddr(
              [&](double x) { return geo(stepper.dense_out(x)).D; }, l, r, 0));
        }
        for (double beta : knots) {
          double ginv = std::sqrt((1 - beta) * (1 + beta));
          auto bound = [beta, ginv](double x, double mu) {
            return x * ginv + beta * mu - 1;
          };
          if (bound(gl.x, gl.mu) * bound(gr.x, gr.mu) <= 0) {
            cuts.push_back(zriddr(
                [&](double x) {
                  auto g = geo(stepper.dense_out(x));
                  return bound(g.x, g.mu);
                },
                l, r, 0));
          }
        }
        l = r;
        gl = gr;
      }
      std::sort(cuts.begin(), cuts.end());
      cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
      for (int i = 0; i < cuts.size() - 1; ++i) {
        double l = cuts[i], r = cuts[i + 1];
        if (geo(stepper.dense_out(std::midpoint(l, r))).D <= 0)
          continue;
        auto integrand = [&](double path_length) {
          auto g = geo(stepper.dense_out(path_length));
          std::array<double, 2> betas;
          if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu,
                               (1 + g.x) * (1 - g.x), betas))
            return 0.0;
          auto rates = rate(betas);
          return (rates[0] + rates[1]) * g.pref;
        };
        double dtau = quad.qags(integrand, l, r, quad_atol, quad_rtol);

        if (tau < 0 && tau + dtau >= 0) {
          double scattered_point = r;
          if (tau + dtau > 0) {
            auto func = [&](double r) {
              if (l == r)
                return tau;
              return tau + quad.qags(integrand, l, r, quad_atol, quad_rtol);
            };
            auto const &dfunc = integrand;
            scattered_point = rtsafe(func, dfunc, l, r, orbit_atol);
          }

          auto r_psi_alpha = stepper.dense_out(scattered_point);
          auto g = geo(r_psi_alpha);
          std::array<double, 2> betas;
          if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu,
                               (1 + g.x) * (1 - g.x), betas)) {
            throw std::runtime_error("No valid beta at scattering point");
          }
          auto rates = rate(betas);
          double sum_rates = rates[0] + rates[1];
          if (sum_rates <= 0.0 || !std::isfinite(sum_rates)) {
            throw std::runtime_error("Invalid scattering rates");
          }

          scattered_data.r_psi_alpha = r_psi_alpha;
          scattered_data.beta =
              (ran.U() * sum_rates < rates[0]) ? betas[0] : betas[1];
          return EvolveResult::Scattered;
        }
        tau += dtau;
      }

      if (event_id == 0) {
        return EvolveResult::Absorbed;
      } else if (event_id == 1) {
        auto [r, psi, alpha] = stepper.y_new;
        double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
        double3 psi_hat = cross(n, r_hat);
        double3 k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
        escaped_data.omega_inf = omega_inf;
        escaped_data.muk = std::clamp(k_hat.z, -1.0, 1.0);
        return EvolveResult::Escaped;
      }
      stepper.update_old();
    }
  }

  void perform_scattering() {
    auto [r, psi, alpha] = scattered_data.r_psi_alpha;
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double3 B_vec = bfield.calc_B(r, std::clamp(r_hat.z, -1.0, 1.0));
    double B = B_vec.length();
    double3 b = B_vec / B;
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 theta_hat = rho > 0 ? double3{r_hat.x * r_hat.z / rho,
                                          r_hat.y * r_hat.z / rho, -rho}
                                : double3{std::copysign(1.0, r_hat.z), 0, 0};
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu = std::clamp(
        b.x * std::cos(alpha) +
            std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat)),
        -1.0, 1.0);

    double beta = scattered_data.beta;
    double mup_out = sample_mup(ran);
    double mu_out =
        std::clamp((mup_out + beta) / (1 + beta * mup_out), -1.0, 1.0);
    double3 b_cartesian = b.x * r_hat + b.y * theta_hat + b.z * phi_hat;
    double3 t_hat = ran.unit_perp_to(b_cartesian);
    double3 k_out =
        mu_out * b_cartesian + std::sqrt(1 - mu_out * mu_out) * t_hat;
    double cross_len = cross(r_hat, k_out).length();
    double alpha_out = std::atan2(cross_len, dot(k_out, r_hat));
    double3 n_out = (cross_len < 1e-10) ? ran.unit_perp_to(r_hat)
                                        : cross(r_hat, k_out) / cross_len;
    double3 e1_out = r_hat;
    double3 e2_out = cross(n_out, e1_out);
    double omega_inf_out = omega_inf * (1 - beta * mu) / (1 - beta * mu_out);
    Polarization pol_out = (ran.U() < 1 / (1 + mup_out * mup_out))
                               ? Polarization::E
                               : Polarization::O;

    if (!std::isfinite(omega_inf_out)) {
      throw std::runtime_error(
          "Non-finite omega_inf encountered after scattering");
    }

    init(n_out, e1_out, e2_out, omega_inf_out, pol_out, r, 0, alpha_out);
  }
};

BField bfield("table/bfield_t10.txt", B_pole, R_star);
Boltzmann fb(-0.75);

int main() {
  PhotonEvolution pe(bfield, fb, 7774);
  constexpr int N = 1e4;
  FILE *f = std::fopen("output/main.txt", "w");
  for (int i = 0; i < N; ++i) {
    pe.init_random_radial(1.0, Polarization::E);
    while (true) {
      auto result = pe.evolve_geodesic();
      if (result == PhotonEvolution::EvolveResult::Escaped) {
        auto [omega_inf, muk] = pe.escaped_data;
        // printf("Escaped: omega_inf = %g, muk = %g\n", omega_inf, muk);
        std::println(f, "{} {}", omega_inf, muk);
        std::print("{}/{}\r", i + 1, N);
        break;
      } else if (result == PhotonEvolution::EvolveResult::Absorbed) {
        // printf("Photon absorbed by the star.\n");
        break;
      } else if (result == PhotonEvolution::EvolveResult::Scattered) {
        pe.perform_scattering();
        // printf("Photon scattered. New state initialized.\n");
      }
    }
  }
  std::fclose(f);
}
