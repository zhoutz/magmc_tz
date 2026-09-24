// Independent reference from doc/bench_optical_depth.typ: integrate in beta,
// explicitly split the admissible domain BEFORE asking QAGS to integrate it.
#include "bench_driver.hpp"

double velocity_truth(BField const &field,transport::Distribution const &dist,Photon const &photon) {
  using namespace transport;
  Medium medium{field,dist,photon};auto point=medium.at({photon.r,0,0});
  double mu=point.mu, q=2+field.p;
  auto x=[&](double r){return B_to_omega*field.calc_B(r,photon.e1.z).length()*std::sqrt(1-rs/r)/photon.omega_inf;};
  double xlo=x(10000),xhi=x(photon.r);
  auto integrand=[&](double beta) {
    double g=(1-beta*mu)/std::sqrt((1-beta)*(1+beta));
    if(!(g>xlo && g<xhi))return 0.;
    double r=zriddr([&](double r){return x(r)-g;},photon.r,10000,1e-11);
    double lapse=std::sqrt(1-rs/r),Q=q-.5*rs/(r-rs);
    double overlap=photon.pol==Polarization::E ? .5 : .5*std::pow((mu-beta)/(1-beta*mu),2);
    return dist.f(beta)*(1-beta*mu)*overlap/(lapse*Q);
  };
  auto knots=dist.knots;knots.push_back(dist.low);knots.push_back(dist.high);
  // Roots at the radial endpoints delimit all admissible velocity intervals,
  // including both sides of the minimum g(beta=mu) in the southern hemisphere.
  for(double edge:{xlo,xhi}) {
    ResonancePoint p{std::log(edge),mu,0};
    for(double beta:medium.branches(p,p.discriminant()))
      if(beta>dist.low && beta<dist.high)knots.push_back(beta);
  }
  std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
  double integral=0;
  for(size_t i=1;i<knots.size();++i)
    integral+=qags(integrand,knots[i-1],knots[i],1e-12,1e-11);
  return (field.p+1)*pi*field.Bphi_over_Btheta(photon.e1.z)/std::abs(dist.mean)*integral;
}

int main() {
  using namespace transport;
  BField field("table/bfield_t10.txt",B_pole,R_star);
  std::map<double,Distribution> distributions;
  auto cases=benchmark_cases("radial");
  for(auto c:cases)if(!distributions.count(c.b0))distributions.emplace(c.b0,boltzmann_distribution(c.b0));
  std::ofstream out("output/radial_truth.txt");out<<std::setprecision(17);
  for(auto c:cases) {
    auto p=make_photon(10,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
    out<<c.b0<<' '<<c.muz<<' '<<c.energy<<' '<<c.pol<<' '<<velocity_truth(field,distributions.at(c.b0),p)<<'\n';
  }
}
