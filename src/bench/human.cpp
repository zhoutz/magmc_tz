#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
// #include "../qags.hpp"
#include "../quad.hpp"
#include "../solve_quadratic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <print>

constexpr double M_star = 1.4;                                    // M_sun
constexpr double R_star = 10;                                     // km
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km; // km
constexpr double B_pole = 1e14;                                   // G

using YVector = std::array<double, 3>;

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  Quad quad;
  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  void operator()(double x, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha] = y;
    double L = std::sqrt(1 - rs / r);
    dydx[0] = L * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * L) * (1 - 3 * rs / (2 * r));
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
    double3 theta_hat = rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double mu = std::clamp(b.x * std::cos(alpha) +
                               std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat)),
                           -1.0, 1.0);
    return Geometry{
        .x = x,
        .mu = mu,
        .D = std::fma(x, x, (mu - 1) * (mu + 1)),
        .pref = (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz) / (std::abs(fb.b_mean) * r),
    };
  }
};

double event_absorb(double x, YVector const &y) { return y[0] - R_star; }
double event_escape(double x, YVector const &y) { return y[0] - 10000; }

BField bfield("table/bfield_t10.txt", B_pole, R_star);

double total_optical_depth(double b0, double muz, double oi, Polarization pol, int n_knots = 8,
                           int n_scan = 8, double atol = 1e-10, double rtol = 1e-9) {
  Boltzmann fb(b0, n_knots);
  double3 r_hat{std::sqrt(1 - muz * muz), 0, muz};
  double3 n{0, 1, 0};

  PhotonEvolution pe{
      .bfield = bfield,
      .fb = fb,
      .n = n,
      .e1 = r_hat,
      .e2 = cross(n, r_hat),
      .omega_inf = oi,
      .pol = pol,
  };

  StepperDopr5<3, PhotonEvolution> stepper(pe, atol, rtol);
  stepper.add_event(&event_absorb);
  stepper.add_event(&event_escape);

  stepper.init(0.0, 1e-3 * R_star, {R_star, 0.0, 0.0});
  double tau = 0;

  while (true) {
    stepper.do_step(0.1 * stepper.y_old[0]);
    int event_id = stepper.detect_event();

    {
      double xl = stepper.x_old, xr = stepper.x_old + stepper.h_old;
      std::vector<double> cuts{xl, xr};
      double l = xl;
      auto gl = pe.geo(stepper.dense_out(l));
      for (int i = 1; i <= n_scan; ++i) {
        double r = xl + (xr - xl) * i / n_scan;
        auto gr = pe.geo(stepper.dense_out(r));
        if (gl.D * gr.D <= 0) {
          cuts.push_back(zriddr([&](double x) { return pe.geo(stepper.dense_out(x)).D; }, l, r, 0));
        }
        for (double beta : fb.knots) {
          auto bound = [beta](double x, double mu) {
            return x * std::sqrt((1 - beta) * (1 + beta)) + beta * mu - 1;
          };
          if (bound(gl.x, gl.mu) * bound(gr.x, gr.mu) <= 0) {
            cuts.push_back(zriddr(
                [&](double x) {
                  auto g = pe.geo(stepper.dense_out(x));
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
        double l = cuts[i], r = cuts[i + 1], m = std::midpoint(l, r);
        if (pe.geo(stepper.dense_out(m)).D <= 0) continue;
        auto integrand = [&](double path_length) {
          auto y = stepper.dense_out(path_length);
          auto g = pe.geo(y);
          std::array<double, 2> betas;
          if (!solve_quadratic(g.x * g.x + g.mu * g.mu, -2 * g.mu, (1 + g.x) * (1 - g.x), betas))
            return 0.0;
          double ret = 0;
          for (double beta : betas) {
            double f = fb.f(beta);
            if (f == 0) continue;
            double overlap = (pol == Polarization::E) ? (0.5) : (0.5 * g.D / (g.x * g.x));
            ret += f * overlap * (1 - beta * g.mu) * (1 + beta) * (1 - beta) * g.x / std::sqrt(g.D);
          }
          return ret * g.pref;
        };
        tau += pe.quad.qags(integrand, l, r, atol, rtol);
      }
    }

    if (event_id != -1) break;
    stepper.update_old();
  }
  return tau;
};

int main() {
  std::FILE *fp = std::fopen("output/human.txt", "w");
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
