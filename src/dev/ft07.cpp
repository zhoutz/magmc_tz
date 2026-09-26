// FT07 Fig. 4: standalone driver derived from main.cpp; existing sources are
// unchanged. Build/run from the repository root; see --help. Energies (omega)
// are in keV. g++-16 src/dev/ft07.cpp -o build/ft07 -std=c++23 -O3 -lgsl
//     -I/opt/homebrew/include -L/opt/homebrew/lib
// build/ft07 --photons 100000
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
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <print>
#include <string>

constexpr double R_star = 10;   // km
constexpr double rs = 0;        // FT07 neglects gravity (sections 3.1 and 3.3).
constexpr double B_pole = 1e14; // G
constexpr double beta0 = 0.75;
constexpr double kT_bb =
    0.4; // keV; chosen physical scale, not a fitted temperature.
constexpr double omega_in = 2.82144 * kT_bb; // FT07 eq. (49).
constexpr double escape_radius =
    10000; // km; convergence can be checked separately.

// FT07 section 3.4 removes 0.001 of the number distribution at each end.
// Keep the full-distribution current normalization: the truncation is a
// numerical approximation, not a new, renormalized particle population.
struct FT07Boltzmann {
  Boltzmann full;
  double b_min, b_max, inv_abs_b_mean;
  FT07Boltzmann(double b0) : full(b0), inv_abs_b_mean(full.inv_abs_b_mean) {
    Quad q;
    auto density = [&](double b) { return full.f(b); };
    auto quantile = [&](double probability) {
      auto residual = [&](double b) {
        if (b == 0)
          return -probability;
        return q.qags(density, 0, b, 1e-11, 1e-10) - probability;
      };
      return zriddr(residual, 0, 1, 1e-12);
    };
    b_min = quantile(0.001);
    b_max = quantile(0.999);
  }
  double f(double b) const { return (b > b_min && b < b_max) ? full.f(b) : 0; }
};

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
  static double event_escape(double, YVector const &y) {
    return y[0] - escape_radius;
  }
  static auto compute_knots(FT07Boltzmann const &fb, int n_knots) {
    std::vector<double> knots(n_knots);
    for (int i = 0; i < n_knots; ++i) {
      knots[i] = fb.b_min + (fb.b_max - fb.b_min) * i / (n_knots - 1);
    }
    return knots;
  }

  BField const &bfield;
  FT07Boltzmann const &fb;
  StepperDopr5<3, decltype(geodesic)> stepper;
  Ran ran;

  std::vector<double> knots;
  const int n_scan;
  const double orbit_atol, orbit_rtol, quad_atol, quad_rtol;

  Quad quad;

  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  PhotonEvolution(BField const &bfield, FT07Boltzmann const &fb, int seed, //
                  int n_knots = 4, int n_scan = 4,                         //
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

  enum class EvolveResult { Surface, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
  } escaped_data;

  struct ScatteredData {
    YVector r_psi_alpha;
    double beta;
  } scattered_data;

  EvolveResult evolve_geodesic() {
    double uniform;
    do {
      uniform = ran.U();
    } while (uniform <= 0 || uniform >= 1);
    double tau = std::log(uniform);
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
        return EvolveResult::Surface;
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

  // FT07 section 3.3: E is absorbed; O has albedo 1/2.
  // The paper does not specify a full surface angular redistribution law.
  // We use the low-frequency magnetized Thomson O->O dipole law,
  // dP/dOmega proportional to 1-(k.B)^2, conditioned on outward escape.
  // Reflection is elastic, and preserves O polarization. This is an explicit
  // surface-model assumption; it is separate from resonant scattering below.
  bool reflect_surface() {
    if (pol == Polarization::E || ran.U() >= 0.5)
      return false;
    auto [r, psi, alpha] = stepper.y_new;
    double3 r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double muz = std::clamp(r_hat.z, -1.0, 1.0);
    double rho = std::hypot(r_hat.x, r_hat.y);
    double3 th = rho > 0
                     ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                     : double3{std::copysign(1.0, muz), 0, 0};
    double3 ph =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    double3 b_sph = to_unit(bfield.calc_B(R_star, muz));
    double3 b = b_sph.x * r_hat + b_sph.y * th + b_sph.z * ph;
    double3 k;
    do {
      k = ran.point_on_unit_sphere();
      if (dot(k, r_hat) <= 0)
        continue;
      double mu = std::clamp(dot(k, b), -1.0, 1.0);
      if (ran.U() < (1 - mu) * (1 + mu))
        break;
    } while (true);
    double3 normal = cross(r_hat, k);
    double sine = normal.length();
    normal = sine > 1e-12 ? normal / sine : ran.unit_perp_to(r_hat);
    // Exact surface position gives the stepper a zero event sign and lets the
    // outward ray depart, instead of detecting the same hit a second time.
    init(normal, r_hat, cross(normal, r_hat), omega_inf, Polarization::O,
         R_star, 0, std::atan2(sine, dot(k, r_hat)));
    return true;
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

int main(int argc, char **argv) {
  try {
    long long N = 100000;
    int seed = 7774;
    std::string output = "output/ft07_beta075_twist10_E.txt";
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--help") {
        std::println(
            "Usage: build/ft07 [--photons N] [--seed SEED] [--output PATH]");
        std::println("Run from the repository root. Defaults: 100000 photons, "
                     "seed 7774.");
        return 0;
      }
      if (i + 1 >= argc)
        throw std::runtime_error("Missing value for " + arg);
      std::string value = argv[++i];
      if (arg == "--photons" || arg == "--seed") {
        std::size_t end;
        if (arg == "--photons")
          N = std::stoll(value, &end);
        else
          seed = std::stoi(value, &end);
        if (end != value.size())
          throw std::runtime_error(arg + " requires a decimal integer");
      } else if (arg == "--output")
        output = value;
      else
        throw std::runtime_error("Unknown option " + arg);
    }
    if (N <= 0)
      throw std::runtime_error("--photons must be positive");
    BField bfield("table/bfield_t10.txt", B_pole, R_star);
    if (std::abs(bfield.Delta_phi - 1) > 1e-10)
      throw std::runtime_error("Expected a Delta_phi = 1 magnetic field table");
    FT07Boltzmann fb(beta0);
    PhotonEvolution pe(bfield, fb, seed);
    auto parent = std::filesystem::path(output).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    std::ofstream out(output);
    if (!out)
      throw std::runtime_error("Cannot open " + output);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << std::setprecision(17)
        << "# FT07 monoenergetic response; only escaped photons are listed\n"
        << "# n_emitted = " << N << '\n'
        << "# omega_in_keV = " << omega_in << '\n'
        << "# kT_bb_keV = " << kT_bb << '\n'
        << "# beta0 = " << beta0 << '\n'
        << "# beta_min = " << fb.b_min << '\n'
        << "# beta_max = " << fb.b_max << '\n'
        << "# beta_mean_full = " << fb.full.b_mean << '\n'
        << "# velocity_retained_probability = 0.998\n"
        << "# delta_phi = " << bfield.Delta_phi << '\n'
        << "# field_p = " << bfield.p << '\n'
        << "# field_table = table/bfield_t10.txt\n"
        << "# seed = " << seed << '\n'
        << "# seed_pol = E\n# seed_direction = radial\n"
        << "# rs_km = " << rs << '\n'
        << "# B_pole_G = " << B_pole << '\n'
        << "# R_star_km = " << R_star << '\n'
        << "# escape_radius_km = " << escape_radius << '\n'
        << "# surface_O_albedo = 0.5\n"
        << "# surface_O_law = outward_sin2_kB_elastic_O_to_O\n"
        << "# orbit_atol = " << pe.orbit_atol << '\n'
        << "# orbit_rtol = " << pe.orbit_rtol << '\n'
        << "# quad_atol = " << pe.quad_atol << '\n'
        << "# quad_rtol = " << pe.quad_rtol << '\n'
        << "# n_scan = " << pe.n_scan << '\n'
        << "# columns: omega_out_keV muk n_scatter n_surface_reflect pol_E\n";
    long long escaped = 0, absorbed = 0, unscattered = 0, surface_hits = 0;
    long long total_reflections = 0, total_scatterings = 0;
    auto start = std::chrono::steady_clock::now();
    for (long long i = 0; i < N; ++i) {
      pe.init_random_radial(omega_in, Polarization::E);
      int n_scatter = 0, n_reflect = 0;
      while (true) {
        auto result = pe.evolve_geodesic();
        if (result == PhotonEvolution::EvolveResult::Escaped) {
          auto [omega, muk] = pe.escaped_data;
          if (!(omega > 0 && std::isfinite(omega) && std::isfinite(muk)))
            throw std::runtime_error("Invalid escaped photon");
          out << omega << ' ' << muk << ' ' << n_scatter << ' ' << n_reflect
              << ' ' << (pe.pol == Polarization::E) << '\n';
          ++escaped;
          if (n_scatter == 0)
            ++unscattered;
          break;
        }
        if (result == PhotonEvolution::EvolveResult::Surface) {
          ++surface_hits;
          if (!pe.reflect_surface()) {
            ++absorbed;
            break;
          }
          ++n_reflect;
          ++total_reflections;
        } else {
          pe.perform_scattering();
          ++n_scatter;
          ++total_scatterings;
        }
      }
      if ((i + 1) % 1000 == 0 || i + 1 == N) {
        double seconds = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();
        std::println(stderr, "{}/{}; escaped={}, absorbed={}; {:.1f} s", i + 1,
                     N, escaped, absorbed, seconds);
        out.flush();
      }
    }
    double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    out << "# n_escaped = " << escaped << '\n'
        << "# n_absorbed = " << absorbed << '\n'
        << "# n_unscattered_escaped = " << unscattered << '\n'
        << "# n_surface_hits = " << surface_hits << '\n'
        << "# n_surface_reflections = " << total_reflections << '\n'
        << "# n_scatterings_total = " << total_scatterings << '\n'
        << "# elapsed_seconds = " << seconds << '\n'
        << "# complete = 1\n";
    out.close();
    std::println("Saved {} ({} escaped / {} emitted)", output, escaped, N);
  } catch (std::exception const &e) {
    std::println(stderr, "FT07 failed: {}", e.what());
    return 1;
  }
  return 0;
}
