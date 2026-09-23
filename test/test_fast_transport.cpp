#include "benchmark_cases.hpp"
#include "../ode/resonance_transport.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>

void require(bool good, std::string const &message) {
  if (!good) throw std::runtime_error(message);
}

int main() {
  try {
    BField field("table/bfield_t10.txt", B_pole, R_star);
    Ran seeds(938475610);
    double max_error = 0, max_position_error = 0, max_tail_difference = 0;
    unsigned position_checks = 0;
    for (int i=0; i<96; ++i) {
      double const temperatures[]{-.75,-.2,-.05,-.005,.2};
      Boltzmann fb(temperatures[i%5]);
      auto p = make_photon(seeds.U(-.999,.999), std::acos(seeds.U(-1.,1.)),
                          seeds.U(0,2*pi), std::exp(seeds.U(std::log(.1),std::log(100.))),
                          i%2 ? Polarization::O : Polarization::E,
                          std::exp(seeds.U(std::log(R_star),std::log(250.))));
      YVector y{p.r,p.psi,p.alpha,0};
      PhotonEvolution physics{field,fb,seeds,p.n,p.e1,p.e2,y,p.omega_inf,p.pol};
      TransportOptions reference_options;
      reference_options.geometry_rtol = 1e-11;
      reference_options.geometry_atol = 1e-12;
      reference_options.quadrature_rtol = 1e-9;
      reference_options.quadrature_atol = 1e-10;
      auto const reference = transport(physics,y,reference_options);
      auto const fast = transport_fast(physics,y);
      double error = std::abs(fast.tau-reference.tau)/std::max(1e-8,reference.tau);
      max_error = std::max(max_error,error);
      require(error<2e-6,"random ray "+std::to_string(i)+" depth error="+std::to_string(error));
      require(fast.termination==reference.termination,"random ray boundary mismatch");
      if (i%8==0) {
        TransportOptions no_tail;
        no_tail.fast_tail_u = std::numeric_limits<double>::infinity();
        no_tail.fast_split_field_knots = true;
        auto const check = transport_fast(physics,y,no_tail);
        double difference = std::abs(check.tau-fast.tau)/std::max(1e-8,reference.tau);
        max_tail_difference = std::max(max_tail_difference,difference);
        require(difference<2e-6,"tail / field-knot convergence failed");
      }
      if (reference.tau>1e-6 && position_checks<24) {
        for(double fraction : {.1,.5,.9}) {
          auto old_stop = transport(physics,y,reference_options,fraction*reference.tau);
          auto new_stop = transport_fast(physics,y,{},fraction*reference.tau);
          require(old_stop.termination==TransportTermination::scattered &&
                  new_stop.termination==TransportTermination::scattered,"lost random scattering event");
          double distance_error = std::abs(old_stop.distance-new_stop.distance)/
                                   std::max(1.,old_stop.distance);
          max_position_error = std::max(max_position_error,distance_error);
          require(distance_error<2e-6,"random scattering position error");
          require(physics.rate(physics.geometry(new_stop.state))>0,"scattering outside supported layer");
        }
        ++position_checks;
      }
    }
    std::cout << std::setprecision(9) << "Fast held-out: 96 random rays, " << position_checks*3
              << " scattering targets; max depth error=" << max_error
              << ", max distance error=" << max_position_error
              << ", tail/knot convergence=" << max_tail_difference << '\n';
  } catch(std::exception const &error) {
    std::cerr << "Fast transport regression failed: " << error.what() << '\n';
    return 1;
  }
}
