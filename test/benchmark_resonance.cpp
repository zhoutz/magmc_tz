#include "../ode/transport_physics.hpp"
#include "../ode/resonance_limiter.hpp"
#include "../ode/resonance_transport.hpp"
#include "../ode/dopr5.hpp"
#include "resonance_reference.hpp"
#include "benchmark_cases.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct Measurement {
  double tau = 0;
  std::size_t steps = 0, rhs = 0, field = 0, rate = 0, rejected = 0, roots = 0;
  std::string status = "ok";
};

Measurement rk(PhotonEvolution const &physics, YVector y, std::string const &method) {
  double tol = method == "rk" || method == "cap_1pct" ? 1e-6 : 1e-8;
  StepperDopr5<4, PhotonEvolution> stepper(physics, tol, tol);
  ThermalStepLimiter thermal{physics};
  if (method == "thermal") stepper.set_step_limiter(thermal);
  else if (method != "rk") stepper.set_step_limiter([&](double, YVector const &v) {
    return (method == "cap_1pct" ? .01 : .001) * v[0];
  });
  stepper.add_event([](double, YVector const &v) { return v[0] - 1000 * R_star; });
  stepper.add_event([](double, YVector const &v) { return v[0] - R_star; });
  // Same order-one residual optical-depth scale as exponential MC thresholds.
  y[3] = -1;
  stepper.init(0, .01, y);
  Measurement ret;
  try {
    for (std::size_t k = 0; k < 1000000; ++k) {
      if (std::cos(stepper.y_old[2]) < 0)
        stepper.h_old = std::min(stepper.h_old, .25 * (stepper.y_old[0] - rs));
      stepper.do_step();
      int event = stepper.detect_event();
      if (event >= 0) { ret.tau = stepper.y_new[3] + 1; break; }
      if (k == 999999) throw std::runtime_error("maximum steps");
      stepper.update_old();
    }
  } catch (std::exception const &e) { ret.status = e.what(); ret.tau = std::numeric_limits<double>::quiet_NaN(); }
  ret.steps = stepper.stats.accepted_steps;
  ret.rhs = stepper.stats.derivative_calls;
  ret.rejected = stepper.stats.rejected_steps;
  ret.rate = ret.rhs;
  ret.field = ret.rhs + (method == "thermal" ? 2 * stepper.stats.limiter_calls : 0);
  return ret;
}

Measurement measure(Case const &c, BField const &field, std::string const &method) {
  Boltzmann fb(c.b0);
  Ran random(1234);
  auto const &p = c.photon;
  YVector y{p.r, p.psi, p.alpha, 0};
  PhotonEvolution physics{field, fb, random, p.n, p.e1, p.e2, y, p.omega_inf, p.pol};
  if (method == "rk" || method.starts_with("cap") || method == "thermal") return rk(physics, y, method);
  TransportOptions opt;
  if (method == "event_tight" || method == "reference" || method == "reference_finer") {
    opt.geometry_rtol = 1e-11;
    opt.geometry_atol = 1e-12;
    opt.quadrature_rtol = 1e-9;
    opt.quadrature_atol = 1e-10;
  }
  if (method == "reference" || method == "reference_finer") {
    opt.max_step_fraction = method == "reference" ? .025 : .0125;
    opt.geometry_rtol = 2e-13;
    opt.geometry_atol = 2e-14;
    opt.quadrature_rtol = 2e-10;
    opt.quadrature_atol = 2e-11;
    if (method == "reference_finer") opt.discriminant_fit_step = 1.25e-4;
  }
  Measurement ret;
  try {
    auto result = method == "fast" ? transport_fast(physics, y, opt) : transport(physics, y, opt);
    ret.tau = result.tau;
    ret.steps = result.stats.geometry_steps;
    ret.rhs = result.stats.geometry_rhs_evaluations;
    ret.rate = result.stats.rate_evaluations;
    ret.roots = result.stats.event_roots;
    ret.field = result.stats.field_evaluations;
  } catch (std::exception const &e) { ret.status = e.what(); ret.tau = std::numeric_limits<double>::quiet_NaN(); }
  return ret;
}

int main(int argc, char **argv) {
  bool quick = false, fast_comparison = false;
  int repeats = 5;
  std::string output = "benchmark/results.csv";
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--quick") quick = true;
    else if (arg == "--fast-comparison") fast_comparison = true;
    else if (arg == "--repeats" && i + 1 < argc) repeats = std::stoi(argv[++i]);
    else if (arg == "--output" && i + 1 < argc) output = argv[++i];
    else throw std::runtime_error("usage: benchmark_resonance [--quick] [--fast-comparison] [--repeats N] [--output path]");
  }
  if (repeats < 1) throw std::runtime_error("repeats must be positive");
  BField field("table/bfield_t10.txt", B_pole, R_star);
  std::ofstream out(output);
  if (!out) throw std::runtime_error("cannot open benchmark output");
  out << std::setprecision(17);
  out << "case,radial,beta0,muz,energy_keV,r0,alpha,azimuth,pol,method,tau_reference,reference_uncertainty,tau,absolute_error,relative_error,probability_error,median_us,steps,rhs_calls,field_calls,rate_calls,rejected,event_roots,status\n";
  auto suite = cases(quick);
  int index = 0;
  for (auto const &c : suite) {
    double reference, uncertainty;
    if (c.radial) {
      auto r = resonance_reference::radial(field, Boltzmann(c.b0), rs, c.photon,
                  std::cos(c.photon.alpha) > 0 ? 1000 * R_star : R_star, 1e-11);
      reference = r.tau;
      auto r2 = resonance_reference::radial(field, Boltzmann(c.b0), rs, c.photon,
                  std::cos(c.photon.alpha) > 0 ? 1000 * R_star : R_star, 1e-12);
      uncertainty = std::max(r.error, std::abs(r.tau - r2.tau));
    } else {
      auto r = measure(c, field, "reference"), r2 = measure(c, field, "reference_finer");
      if (r.status != "ok" || r2.status != "ok") {
        std::cerr << "reference failed " << c.name << ": " << r.status << " / " << r2.status << '\n';
        reference = std::numeric_limits<double>::quiet_NaN();
        uncertainty = reference;
      } else { reference = r2.tau; uncertainty = std::abs(r.tau - r2.tau); }
    }
    auto methods = fast_comparison ? std::vector<std::string>{"event", "fast"} :
        std::vector<std::string>{"rk", "cap_1pct", "cap_0.1pct", "thermal", "event", "event_tight", "fast"};
    if (fast_comparison && index%2) std::reverse(methods.begin(), methods.end());
    for (std::string method : methods) {
      std::vector<double> times;
      auto m = measure(c, field, method); // untimed warmup
      for (int j = 0; j < repeats; ++j) {
        auto begin = std::chrono::steady_clock::now();
        m = measure(c, field, method);
        times.push_back(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now()-begin).count());
      }
      std::sort(times.begin(), times.end());
      double absolute = std::abs(m.tau - reference);
      double relative = reference > 1e-14 ? absolute / reference : absolute;
      auto const &p = c.photon;
      double az = std::atan2(p.e2.y, p.e2.x * p.e1.z - p.e2.z * p.e1.x);
      out << c.name << ',' << c.radial << ',' << c.b0 << ',' << p.e1.z << ',' << p.omega_inf << ','
          << p.r << ',' << p.alpha << ',' << az << ',' << (p.pol == Polarization::O ? 'O' : 'E') << ',' << method << ','
          << reference << ',' << uncertainty << ',' << m.tau << ',' << absolute << ',' << relative << ','
          << std::abs(std::exp(-m.tau)-std::exp(-reference)) << ',' << times[times.size()/2] << ','
          << m.steps << ',' << m.rhs << ',' << m.field << ',' << m.rate << ',' << m.rejected << ',' << m.roots
          << ',' << '"' << m.status << '"' << '\n';
      out.flush();
    }
    std::cerr << ++index << '/' << suite.size() << ' ' << c.name << " reference=" << reference << '\n';
  }
}
