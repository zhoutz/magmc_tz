#include "bfield.hpp"
#include "constants.hpp"
#include "distribution.hpp"
#include "dopr5.hpp"
#include "dopr5_maxstep.hpp"
#include "dopr5_dtau_limit.hpp"
#include "init.hpp"
#include "photon.hpp"
#include "ran.hpp"
#include "sample_mup.hpp"
#include "solve_quadratic.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <print>
#include <string>
#include <vector>

constexpr double M_star = 1.4;
constexpr double R_star = 10;
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km;
constexpr double B_pole = 1e14;

using YVector = std::array<double, 4>;

struct PhotonEvolution {
  BField const &bfield;
  Boltzmann const &fb;
  Ran &ran;
  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  void operator()(double x, YVector const &y, YVector &dydx) const {
    auto [r, psi, alpha, tau] = y;
    double f = std::sqrt(1 - rs / r);

    dydx[0] = f * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
    dydx[3] = calc_dtaudl(n, e1, e2, r, psi, alpha, omega_inf, pol);
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
    double rho = std::sqrt(r_hat.x * r_hat.x + r_hat.y * r_hat.y);
    double3 theta_hat{r_hat.x * r_hat.z / rho, r_hat.y * r_hat.z / rho, -rho};
    double3 phi_hat{-r_hat.y / rho, r_hat.x / rho, 0};
    double mu_in =
        b.x * std::cos(alpha) + std::sin(alpha) * (b.y * dot(n, phi_hat) - b.z * dot(n, theta_hat));
    std::array<double, 2> betas;
    if (!solve_quadratic(x * x + mu_in * mu_in, -2 * mu_in, 1 - x * x, betas)) return 0;
    double ret = 0;
    for (double beta : betas) {
      double f = fb.f(beta);
      if (f == 0) continue;
      double mup_in = (mu_in - beta) / (1 - beta * mu_in);
      double esq = (pol == Polarization::E) ? (0.5) : (0.5 * mup_in * mup_in);
      ret += f * esq * (1 - beta * mu_in) * (1 - beta * mu_in) * (1 - beta * beta) /
             std::abs(mu_in - beta);
    }
    ret *= (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz) / (std::abs(fb.b_bar()) * r);

    if (!std::isfinite(ret)) {
      throw std::runtime_error("Non-finite dtaudl encountered");
    }

    return ret;
  }
};

struct TestCase {
  std::string name;
  double r_init;
  double alpha_init;
  double omega_inf;
  Polarization pol;
  int seed;
};

struct BenchmarkResult {
  std::string method_name;
  std::string test_case;
  double final_tau;
  int num_steps;
  double time_ms;
  bool escaped;
  bool absorbed;
};

// Method 1: Original (baseline)
BenchmarkResult run_method_original(BField const &bfield, Boltzmann const &fb, TestCase const &tc) {
  auto start = std::chrono::high_resolution_clock::now();

  Ran ran(tc.seed);
  double3 n = {0, 0, 1};
  double3 e1 = {1, 0, 0};
  double3 e2 = {0, 1, 0};

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .ran = ran,
      .n = n,
      .e1 = e1,
      .e2 = e2,
      .omega_inf = tc.omega_inf,
      .pol = tc.pol,
  };

  YVector y = {tc.r_init, 0.0, tc.alpha_init, 0.0};
  StepperDopr5<4, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);

  auto event_escape = [](double x, YVector const &y) { return y[0] - 1000 * R_star; };
  auto event_absorption = [](double x, YVector const &y) { return y[0] - R_star; };

  stepper.add_event(event_escape);
  stepper.add_event(event_absorption);
  stepper.init(0.0, 1e-3 * R_star, y);

  int num_steps = 0;
  bool escaped = false;
  bool absorbed = false;
  double final_tau = 0.0;

  while (num_steps < 1000000) {
    stepper.do_step();
    num_steps++;

    int event_id = stepper.detect_event();
    if (event_id == 0) {
      escaped = true;
      final_tau = stepper.y_new[3];
      break;
    } else if (event_id == 1) {
      absorbed = true;
      final_tau = stepper.y_new[3];
      break;
    }
    stepper.update_old();
  }

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

  return BenchmarkResult{
      .method_name = "Original",
      .test_case = tc.name,
      .final_tau = final_tau,
      .num_steps = num_steps,
      .time_ms = time_ms,
      .escaped = escaped,
      .absorbed = absorbed,
  };
}

// Method 2: Max step size
BenchmarkResult run_method_maxstep(BField const &bfield, Boltzmann const &fb, TestCase const &tc,
                                   double h_max) {
  auto start = std::chrono::high_resolution_clock::now();

  Ran ran(tc.seed);
  double3 n = {0, 0, 1};
  double3 e1 = {1, 0, 0};
  double3 e2 = {0, 1, 0};

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .ran = ran,
      .n = n,
      .e1 = e1,
      .e2 = e2,
      .omega_inf = tc.omega_inf,
      .pol = tc.pol,
  };

  YVector y = {tc.r_init, 0.0, tc.alpha_init, 0.0};
  StepperDopr5MaxStep<4, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6, h_max);

  auto event_escape = [](double x, YVector const &y) { return y[0] - 1000 * R_star; };
  auto event_absorption = [](double x, YVector const &y) { return y[0] - R_star; };

  stepper.add_event(event_escape);
  stepper.add_event(event_absorption);
  stepper.init(0.0, 1e-3 * R_star, y);

  int num_steps = 0;
  bool escaped = false;
  bool absorbed = false;
  double final_tau = 0.0;

  while (num_steps < 1000000) {
    stepper.do_step();
    num_steps++;

    int event_id = stepper.detect_event();
    if (event_id == 0) {
      escaped = true;
      final_tau = stepper.y_new[3];
      break;
    } else if (event_id == 1) {
      absorbed = true;
      final_tau = stepper.y_new[3];
      break;
    }
    stepper.update_old();
  }

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

  return BenchmarkResult{
      .method_name = "MaxStep_" + std::to_string(h_max),
      .test_case = tc.name,
      .final_tau = final_tau,
      .num_steps = num_steps,
      .time_ms = time_ms,
      .escaped = escaped,
      .absorbed = absorbed,
  };
}

// Method 3: dtau/dl limiting
BenchmarkResult run_method_dtau_limit(BField const &bfield, Boltzmann const &fb,
                                      TestCase const &tc, double dtau_max) {
  auto start = std::chrono::high_resolution_clock::now();

  Ran ran(tc.seed);
  double3 n = {0, 0, 1};
  double3 e1 = {1, 0, 0};
  double3 e2 = {0, 1, 0};

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .ran = ran,
      .n = n,
      .e1 = e1,
      .e2 = e2,
      .omega_inf = tc.omega_inf,
      .pol = tc.pol,
  };

  YVector y = {tc.r_init, 0.0, tc.alpha_init, 0.0};
  StepperDopr5DtauLimit<4, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6, dtau_max, 3);

  auto event_escape = [](double x, YVector const &y) { return y[0] - 1000 * R_star; };
  auto event_absorption = [](double x, YVector const &y) { return y[0] - R_star; };

  stepper.add_event(event_escape);
  stepper.add_event(event_absorption);
  stepper.init(0.0, 1e-3 * R_star, y);

  int num_steps = 0;
  bool escaped = false;
  bool absorbed = false;
  double final_tau = 0.0;

  while (num_steps < 1000000) {
    stepper.do_step();
    num_steps++;

    int event_id = stepper.detect_event();
    if (event_id == 0) {
      escaped = true;
      final_tau = stepper.y_new[3];
      break;
    } else if (event_id == 1) {
      absorbed = true;
      final_tau = stepper.y_new[3];
      break;
    }
    stepper.update_old();
  }

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

  return BenchmarkResult{
      .method_name = "DtauLimit_" + std::to_string(dtau_max),
      .test_case = tc.name,
      .final_tau = final_tau,
      .num_steps = num_steps,
      .time_ms = time_ms,
      .escaped = escaped,
      .absorbed = absorbed,
  };
}

// Method 4: Very tight tolerances (reference solution)
BenchmarkResult run_method_reference(BField const &bfield, Boltzmann const &fb,
                                     TestCase const &tc) {
  auto start = std::chrono::high_resolution_clock::now();

  Ran ran(tc.seed);
  double3 n = {0, 0, 1};
  double3 e1 = {1, 0, 0};
  double3 e2 = {0, 1, 0};

  PhotonEvolution photon_evolution{
      .bfield = bfield,
      .fb = fb,
      .ran = ran,
      .n = n,
      .e1 = e1,
      .e2 = e2,
      .omega_inf = tc.omega_inf,
      .pol = tc.pol,
  };

  YVector y = {tc.r_init, 0.0, tc.alpha_init, 0.0};
  StepperDopr5<4, PhotonEvolution> stepper(photon_evolution, 1e-10, 1e-10);

  auto event_escape = [](double x, YVector const &y) { return y[0] - 1000 * R_star; };
  auto event_absorption = [](double x, YVector const &y) { return y[0] - R_star; };

  stepper.add_event(event_escape);
  stepper.add_event(event_absorption);
  stepper.init(0.0, 1e-6 * R_star, y);

  int num_steps = 0;
  bool escaped = false;
  bool absorbed = false;
  double final_tau = 0.0;

  while (num_steps < 10000000) {
    stepper.do_step();
    num_steps++;

    int event_id = stepper.detect_event();
    if (event_id == 0) {
      escaped = true;
      final_tau = stepper.y_new[3];
      break;
    } else if (event_id == 1) {
      absorbed = true;
      final_tau = stepper.y_new[3];
      break;
    }
    stepper.update_old();
  }

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

  return BenchmarkResult{
      .method_name = "Reference",
      .test_case = tc.name,
      .final_tau = final_tau,
      .num_steps = num_steps,
      .time_ms = time_ms,
      .escaped = escaped,
      .absorbed = absorbed,
  };
}

int main() {
  BField bfield("table/bfield_t10.txt", B_pole, R_star);
  Boltzmann fb(-0.75);

  // Define test cases
  std::vector<TestCase> test_cases = {
      {"radial_out_close", 1.5 * R_star, 0.1, 1.0, Polarization::E, 1001},
      {"radial_out_far", 5.0 * R_star, 0.05, 1.0, Polarization::E, 1002},
      {"tangential", 3.0 * R_star, pi / 2 - 0.1, 1.0, Polarization::E, 1003},
      {"oblique_30", 2.5 * R_star, pi / 6, 1.0, Polarization::E, 1004},
      {"oblique_60", 2.5 * R_star, pi / 3, 1.0, Polarization::E, 1005},
      {"high_energy", 2.0 * R_star, 0.2, 5.0, Polarization::E, 1006},
      {"low_energy", 2.0 * R_star, 0.2, 0.5, Polarization::E, 1007},
      {"O_mode", 2.0 * R_star, 0.2, 1.0, Polarization::O, 1008},
  };

  std::vector<BenchmarkResult> all_results;

  std::println("=== COMPREHENSIVE BENCHMARK ===\n");
  std::println("Computing reference solutions (very tight tolerances)...");

  for (auto const &tc : test_cases) {
    std::println("  Reference for {}", tc.name);
    try {
      auto result = run_method_reference(bfield, fb, tc);
      all_results.push_back(result);
      std::println("    tau = {:.8f}, steps = {}, time = {:.1f} ms", result.final_tau,
                   result.num_steps, result.time_ms);
    } catch (std::exception const &e) {
      std::println("    ERROR: {}", e.what());
    }
  }

  std::println("\nTesting Method 1: Original (baseline)...");
  for (auto const &tc : test_cases) {
    std::println("  Original for {}", tc.name);
    try {
      auto result = run_method_original(bfield, fb, tc);
      all_results.push_back(result);
      std::println("    tau = {:.8f}, steps = {}, time = {:.1f} ms", result.final_tau,
                   result.num_steps, result.time_ms);
    } catch (std::exception const &e) {
      std::println("    ERROR: {}", e.what());
    }
  }

  std::println("\nTesting Method 2: Max Step Size...");
  std::vector<double> h_max_values = {0.1 * R_star, 0.05 * R_star, 0.01 * R_star};
  for (double h_max : h_max_values) {
    std::println("  h_max = {} km", h_max);
    for (auto const &tc : test_cases) {
      std::println("    {}", tc.name);
      try {
        auto result = run_method_maxstep(bfield, fb, tc, h_max);
        all_results.push_back(result);
        std::println("      tau = {:.8f}, steps = {}, time = {:.1f} ms", result.final_tau,
                     result.num_steps, result.time_ms);
      } catch (std::exception const &e) {
        std::println("      ERROR: {}", e.what());
      }
    }
  }

  std::println("\nTesting Method 3: Dtau/dl Limiting...");
  std::vector<double> dtau_max_values = {0.1, 0.05, 0.01};
  for (double dtau_max : dtau_max_values) {
    std::println("  dtau_max = {}", dtau_max);
    for (auto const &tc : test_cases) {
      std::println("    {}", tc.name);
      try {
        auto result = run_method_dtau_limit(bfield, fb, tc, dtau_max);
        all_results.push_back(result);
        std::println("      tau = {:.8f}, steps = {}, time = {:.1f} ms", result.final_tau,
                     result.num_steps, result.time_ms);
      } catch (std::exception const &e) {
        std::println("      ERROR: {}", e.what());
      }
    }
  }

  // Save all results
  std::ofstream fout("benchmark_comprehensive.txt");
  fout << std::setprecision(10);
  fout << "method\ttest_case\tfinal_tau\tnum_steps\ttime_ms\tescaped\tabsorbed\n";
  for (auto const &r : all_results) {
    fout << r.method_name << "\t" << r.test_case << "\t" << r.final_tau << "\t" << r.num_steps
         << "\t" << r.time_ms << "\t" << r.escaped << "\t" << r.absorbed << "\n";
  }

  std::println("\n=== RESULTS SAVED ===");
  std::println("File: benchmark_comprehensive.txt");
  std::println("Total runs: {}", all_results.size());

  return 0;
}
