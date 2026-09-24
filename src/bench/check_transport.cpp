// Supplement total depths with common exponential optical-depth draws and
// cumulative-depth checks. Write each measured discrepancy, without assertions
// that merely repeat implementation choices.
#include "bench_driver.hpp"
int main() {
  using namespace transport;
  gsl_set_error_handler_off();
  BField field("table/bfield_t10.txt",B_pole,R_star);
  auto cases=benchmark_cases("angular");
  std::map<double,Distribution> distributions;
  std::ofstream out("output/scattering_checks.txt");out<<std::setprecision(17);
  out<<"# case method U target ref_tau ref_outcome outcome ref_length length abs_length_error cumulative_error\n";
  for(size_t i=0;i<cases.size();i+=17) {
    auto c=cases[i];
    if(!distributions.count(c.b0))distributions.emplace(c.b0,boltzmann_distribution(c.b0));
    auto const &dist=distributions.at(c.b0);
    auto p=make_photon(c.r,c.muz,c.alpha,c.azimuth,c.energy,c.pol ? Polarization::E : Polarization::O);
    for(double U:{.05,.5,.95}) {
      Options o;o.tolerance=1e-8;o.resolution=.005;o.probes=16;o.tau_target=-std::log(U);
      auto ref=integrate(field,dist,p,Method::reference,o);
      for(auto method:{Method::ft07,Method::panels,Method::geometric}) {
        o.tolerance=1e-6;o.resolution=method==Method::ft07 ? .01 : method==Method::panels ? .2 : .001;o.probes=8;
        auto test=integrate(field,dist,p,method,o);
        // Evaluate the reference cumulative depth at the candidate location.
        Options check;check.tolerance=1e-8;check.resolution=.005;check.probes=16;check.path_limit=test.length;
        auto cumulative=integrate(field,dist,p,Method::reference,check);
        out<<i<<' '<<(method==Method::ft07 ? "ft07" : method==Method::panels ? "velocity_panels" : "geometric")
           <<' '<<U<<' '<<o.tau_target<<' '<<ref.tau<<' '<<ref.outcome<<' '<<test.outcome<<' '
           <<ref.length<<' '<<test.length<<' '<<std::abs(ref.length-test.length)<<' '<<std::abs(cumulative.tau-test.tau)<<'\n';
      }
    }
  }
}
