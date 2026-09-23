#include "transport_physics.hpp"
#include "resonance_transport.hpp"
#include "resonance_limiter.hpp"
#include "init.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

// The coupled RK variants remain available for reproducing the benchmark.
TransportResult coupled_flight(PhotonEvolution const &physics, double target,
                               std::string const &method) {
  StepperDopr5<4, PhotonEvolution> stepper(physics, 1e-8, 1e-8);
  if (method == "thermal") stepper.set_step_limiter(ThermalStepLimiter{physics});
  if (method == "cap") stepper.set_step_limiter([](double, YVector const &v) { return .01 * v[0]; });
  stepper.add_event([](double, YVector const &v) { return v[0] - 1000 * R_star; });
  stepper.add_event([](double, YVector const &v) { return v[0] - R_star; });
  stepper.add_event([](double, YVector const &v) { return v[3]; });
  auto y = physics.r_psi_alpha_tau;
  y[3] = -target;
  stepper.init(0, .01, y);
  for (int k = 0; k < 1000000; ++k) {
    if (std::cos(stepper.y_old[2]) < 0)
      stepper.h_old = std::min(stepper.h_old, .25 * (stepper.y_old[0] - rs));
    stepper.do_step();
    int id = stepper.detect_event();
    if (id >= 0) {
      TransportResult result;
      result.state = stepper.y_new;
      result.tau = result.state[3] += target;
      result.distance = stepper.x_old + stepper.h_old;
      result.termination = id == 0 ? TransportTermination::escaped :
                           id == 1 ? TransportTermination::absorbed : TransportTermination::scattered;
      return result;
    }
    stepper.update_old();
  }
  throw std::runtime_error("coupled flight exceeded step limit");
}

int main(int argc, char **argv) {
  try {
    std::string method = "event";
    int photons = 1;
    unsigned long long seed = 1234;
    double beta0 = -.75, energy = 1;
    Polarization mode = Polarization::O;
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--help") {
        std::cout << "Run from project root. Options: --method event|thermal|cap|rk "
                     "--photons N --seed N --beta0 B --energy keV --mode O|E\n";
        return 0;
      }
      if (i + 1 >= argc) throw std::invalid_argument("missing argument value");
      std::string value = argv[++i];
      if (arg == "--method") method = value;
      else if (arg == "--photons") photons = std::stoi(value);
      else if (arg == "--seed") seed = std::stoull(value);
      else if (arg == "--beta0") beta0 = std::stod(value);
      else if (arg == "--energy") energy = std::stod(value);
      else if (arg == "--mode" && (value == "O" || value == "E"))
        mode = value == "O" ? Polarization::O : Polarization::E;
      else throw std::invalid_argument("unknown option or mode: " + arg);
    }
    if (photons < 1 || !(energy > 0) || !std::isfinite(energy) ||
        (method != "event" && method != "thermal" && method != "cap" && method != "rk"))
      throw std::invalid_argument("invalid simulation options");
    BField field("table/bfield_t10.txt", B_pole, R_star);
    Boltzmann distribution(beta0);
    Ran random(seed);
    int escaped = 0, absorbed = 0;
    std::size_t scatterings = 0;
    std::cout << std::setprecision(12);
    for (int i = 0; i < photons; ++i) {
      auto p = init07(random, R_star, mode);
      PhotonEvolution physics{field, distribution, random, p.n, p.e1, p.e2,
                              {p.r, p.psi, p.alpha, 0}, energy, mode};
      bool done = false;
      for (int interactions = 0; interactions < 100000; ++interactions) {
        double target = -std::log(random.U_open());
        auto result = method == "event" ? transport(physics, physics.r_psi_alpha_tau, {}, target)
                                        : coupled_flight(physics, target, method);
        if (photons == 1) {
          std::cout << (result.termination == TransportTermination::scattered ? "scattered" :
                        result.termination == TransportTermination::escaped ? "escaped" : "absorbed")
                    << " r=" << result.state[0] << " psi=" << result.state[1]
                    << " alpha=" << result.state[2] << " tau=" << result.tau << '\n';
        }
        if (result.termination == TransportTermination::escaped) { ++escaped; done = true; break; }
        if (result.termination == TransportTermination::absorbed) { ++absorbed; done = true; break; }
        ++scatterings;
        physics.r_psi_alpha_tau = result.state;
        physics.perform_scattering();
      }
      if (!done) throw std::runtime_error("photon exceeded interaction limit");
    }
    std::cout << "method=" << method << " photons=" << photons << " escaped=" << escaped
              << " absorbed=" << absorbed << " scatterings=" << scatterings << '\n';
  } catch (std::exception const &e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
}
