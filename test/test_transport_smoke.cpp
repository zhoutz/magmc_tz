#include "resonance_transport.hpp"
#include "init.hpp"
#include "resonance_reference.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
bool use_fast = false;
TransportResult tested_transport(PhotonEvolution const &p, YVector y,
                                 TransportOptions const &options = {},
                                 double target = std::numeric_limits<double>::infinity(),
                                 double outer = 1000*R_star) {
  return use_fast ? transport_fast(p,y,options,target,outer) : transport(p,y,options,target,outer);
}
void require(bool condition, std::string const &message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, std::string const &message) {
  require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
          message + ": got " + std::to_string(actual) + ", expected " + std::to_string(expected));
}

Photon radial_photon(double muz, Polarization mode) {
  double const sine = std::sqrt((1 - muz) * (1 + muz));
  return {{0, 1, 0}, {sine, 0, muz}, {muz, 0, -sine}, R_star, 0, 0, 1, mode};
}

void test_scattering_locations(BField const &field, double muz, Polarization mode, double beta0 = -.2) {
  Boltzmann distribution(beta0);
  Ran random(310);
  Photon const photon = radial_photon(muz, mode);
  YVector const initial{photon.r, photon.psi, photon.alpha, 0};
  PhotonEvolution physics{field, distribution, random, photon.n, photon.e1, photon.e2,
                          initial, photon.omega_inf, photon.pol};
  constexpr double outer = 1000 * R_star;
  auto const reference = resonance_reference::radial(field, distribution, rs, photon, outer, 1e-11);
  double const tolerance = 1e-6 * std::max(1.0, reference.tau);
  auto const full = tested_transport(physics, initial, {}, std::numeric_limits<double>::infinity(), outer);
  require(full.termination == TransportTermination::escaped, "radial ray did not escape");
  near(full.tau, reference.tau, tolerance, "full radial tau disagrees with independent reference");
  double last_radius = R_star;
  for (double fraction : {.1, .5, .9}) {
    double const target = fraction * reference.tau;
    auto const stopped = tested_transport(physics, initial, {}, target, outer);
    require(stopped.termination == TransportTermination::scattered,
            "finite optical-depth target failed to produce scattering event");
    near(stopped.tau, target, 1e-14 * std::max(1.0, target), "reported target optical depth");
    require(stopped.state[0] > last_radius && stopped.state[0] < outer,
            "scattering positions are not ordered by target depth");
    last_radius = stopped.state[0];
    // The reference changes variables to carrier velocity; it does not use the
    // production rate, ODE, root solver, or optical-depth quadrature.
    auto const prefix = resonance_reference::radial(field, distribution, rs, photon,
                                                     stopped.state[0], 1e-11);
    near(prefix.tau, target, tolerance, "scattering radius has incorrect reference optical depth");
    require(physics.rate(physics.geometry(stopped.state)) > 0,
            "scattering event located outside supported resonance layer");
    if (fraction == .5) {
      // Check both retained and reset tau offsets when restarting inside a layer.
      auto const continued = tested_transport(physics, stopped.state, {},
                                        std::numeric_limits<double>::infinity(), outer);
      near(continued.tau, reference.tau, 2 * tolerance, "continued accumulated tau lost additivity");
      YVector tail_initial = stopped.state;
      tail_initial[3] = 0;
      auto const tail = tested_transport(physics, tail_initial, {},
                                  std::numeric_limits<double>::infinity(), outer);
      near(prefix.tau + tail.tau, reference.tau, 2 * tolerance,
           "independent prefix plus restarted tail lost additivity");
      auto const first = tested_transport(physics, initial, {},
                                   std::numeric_limits<double>::infinity(), stopped.state[0]);
      require(first.termination == TransportTermination::escaped,
              "finite radial segment failed to reach its outer boundary");
      near(first.tau + tail.tau, full.tau, 2 * tolerance,
           "two production radial segments disagree with full propagation");
    }
  }
}

void test_cold_fold(BField const &field) {
  Boltzmann distribution(-.001);
  Ran random(42);
  auto p = radial_photon(0, Polarization::E);
  PhotonEvolution physics{field, distribution, random, p.n, p.e1, p.e2,
                          {p.r, p.psi, p.alpha, 0}, p.omega_inf, p.pol};
  auto exact = resonance_reference::radial(field, distribution, rs, p, 10000, 1e-11);
  TransportOptions tight;
  tight.geometry_rtol = 1e-11;
  tight.geometry_atol = 1e-12;
  tight.quadrature_rtol = 1e-9;
  tight.quadrature_atol = 1e-10;
  for (double fit_step : {2.5e-4, 1.25e-4}) {
    tight.discriminant_fit_step = fit_step;
    auto value = tested_transport(physics, physics.r_psi_alpha_tau, tight);
    near(value.tau, exact.tau, 1e-7 * exact.tau, "cold E fold / fit-scale convergence");
  }
}

void test_short_ensemble(BField const &field) {
  Boltzmann distribution(-.2);
  Ran random(20260923);
  unsigned escapes = 0, absorptions = 0, scatterings = 0;
  constexpr unsigned photons = 16;
  for (unsigned i = 0; i < photons; ++i) {
    Photon const p = init07(random, R_star, i % 2 ? Polarization::O : Polarization::E);
    PhotonEvolution physics{field, distribution, random, p.n, p.e1, p.e2,
                            {p.r, p.psi, p.alpha, 0}, p.omega_inf, p.pol};
    bool terminated = false;
    for (unsigned segment = 0; segment < 256; ++segment) {
      double const target = -std::log(random.U_open());
      auto const result = tested_transport(physics, physics.r_psi_alpha_tau, {}, target);
      for (double value : result.state)
        require(std::isfinite(value), "photon ensemble produced non-finite transport state");
      require(result.distance >= 0 && std::isfinite(result.distance) && result.tau >= 0 &&
                  std::isfinite(result.tau),
              "photon ensemble produced invalid path length or optical depth");
      if (result.termination == TransportTermination::escaped) {
        near(result.state[0], 1000 * R_star, 1e-7, "ensemble escape radius");
        ++escapes;
        terminated = true;
        break;
      }
      if (result.termination == TransportTermination::absorbed) {
        near(result.state[0], R_star, 1e-8, "ensemble absorption radius");
        ++absorptions;
        terminated = true;
        break;
      }
      ++scatterings;
      physics.r_psi_alpha_tau = result.state;
      physics.perform_scattering();
      require(physics.omega_inf > 0 && std::isfinite(physics.omega_inf),
              "scattering produced invalid photon energy");
      require(physics.r_psi_alpha_tau[3] == 0, "scattering did not reset optical depth");
      near(physics.n.length(), 1, 1e-12, "scattered orbit plane normal norm");
      near(physics.e1.length(), 1, 1e-12, "scattered radial basis norm");
      near(physics.e2.length(), 1, 1e-12, "scattered tangential basis norm");
      near(dot(physics.n, physics.e1), 0, 1e-12, "scattered orbit basis orthogonality");
    }
    require(terminated, "short photon ensemble exceeded 256 scatterings for one photon");
  }
  require(escapes + absorptions == photons && escapes > 0 && scatterings > 0,
          "ensemble did not exercise scattering and boundary termination");
  std::cout << "Ensemble: " << escapes << " escaped, " << absorptions << " absorbed, "
            << scatterings << " scatterings\n";
}
} // namespace

int main() {
  try {
    BField field("table/bfield_t10.txt", B_pole, R_star);
    for (bool fast : {false, true}) {
    use_fast = fast;
    test_scattering_locations(field, -.2, Polarization::O);
    test_scattering_locations(field, -.2, Polarization::E);
    test_scattering_locations(field, 0, Polarization::E);
    test_scattering_locations(field, 0, Polarization::E, -.001);
    test_cold_fold(field);
    test_short_ensemble(field);
    }
    std::cout << "Transport smoke tests passed\n";
    return 0;
  } catch (std::exception const &error) {
    std::cerr << "Transport smoke test failed: " << error.what() << '\n';
    return 1;
  }
}
