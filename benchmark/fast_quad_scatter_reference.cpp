// Frozen, independently refined GL8/16 reference: unaffected by the GK refactor.
#include "../output/quad_refactor_before/src/transport.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
struct Input {int id;double b0,muz,alpha,azimuth,r,psi,energy;int pol;double reference;int outcome;double initial,target;};
int main() {
 try {
 BField field("table/bfield_t10.txt",B_pole,R_star);
 std::ifstream cases_file("output/quad_refactor_cases.txt");std::string line;
 std::map<int,Input> inputs;std::map<int,transport::Result> totals;
 transport::Options ref;ref.method=transport::Method::event_mesh;ref.resolution=.005;ref.tolerance=1e-8;ref.probes=16;
 auto photon=[](Input const &c){auto p=transport::make_photon(c.r,c.muz,c.alpha,c.azimuth,c.energy,Polarization(c.pol));p.psi=c.psi;return p;};
 std::ofstream truth("output/quad_scatter_truth.csv");truth<<std::setprecision(17)<<"id,total_tau,expected_tau,expected_outcome\n";
 while(std::getline(cases_file,line)) {
  if(line.empty()||line[0]=='#')continue;
  std::istringstream row(line);std::string suite;Input c;
  row>>suite>>c.id>>c.b0>>c.muz>>c.alpha>>c.azimuth>>c.r>>c.psi>>c.energy>>c.pol>>c.reference>>c.outcome>>c.initial>>c.target;
  if(suite!="scatter")continue;
  auto total=transport::integrate(field,Boltzmann(c.b0),photon(c),ref);
  inputs[c.id]=c;totals[c.id]=total;
  bool event=c.target<total.tau+c.initial;
  truth<<c.id<<','<<total.tau<<','<<(event?c.target:c.initial+total.tau)<<','<<(event?2:int(total.outcome))<<'\n';
 }
 if(inputs.size()!=36)throw std::runtime_error("Expected 36 scatter inputs");
 std::ifstream locations("output/quad_scatter_locations.txt");
 std::ofstream out("output/quad_scattering_independent.csv");
 out<<std::setprecision(17)<<"variant,rtol,id,outcome,expected_outcome,cumulative,target_or_total,absolute_residual,reported_tau_error,impact_relative_error,status\n";
 std::string variant;double tol,distance,tau;int id,outcome,count=0,failures=0;
 while(locations>>variant>>tol>>id>>distance>>tau>>outcome) {
  auto const &c=inputs.at(id);auto const &total=totals.at(id);bool expected_scatter=c.target<total.tau+c.initial;
  int expected=expected_scatter?2:int(total.outcome);
  double target=expected_scatter?c.target:total.tau+c.initial;
  auto limited=ref;limited.path_limit=distance;
  auto r=transport::integrate(field,Boltzmann(c.b0),photon(c),limited);
  double cumulative=c.initial+r.tau,error=std::abs(cumulative-target);
  double b1=c.r*std::sin(c.alpha)/std::sqrt(1-rs/c.r),b2=r.state[0]*std::sin(r.state[2])/std::sqrt(1-rs/r.state[0]);
  double impact=std::abs(b1-b2)/b1;
  bool ok=outcome==expected && error<=1e-6+1e-5*std::abs(target);
  out<<variant<<','<<tol<<','<<id<<','<<outcome<<','<<expected<<','<<cumulative<<','<<target<<','<<error<<','<<std::abs(tau-target)<<','<<impact<<','<<(ok?"ok":"mismatch")<<'\n';++count;failures+=!ok;
 }
 std::cout<<"Independent scatter-location checks="<<count<<" mismatches="<<failures<<'\n';
 return failures?1:0;
 }catch(std::exception const &e){std::cerr<<e.what()<<'\n';return 1;}
}
