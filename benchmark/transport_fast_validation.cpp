#include "benchmark_common.hpp"
#include "../src/transport_fast.hpp"
#include <random>

namespace {
struct AngularCase {int id;double b0,muz,azimuth;Photon photon;};
struct Answer {double tau,distance;int outcome;long fields,steps;std::string status="ok";};
using Clock=std::chrono::steady_clock;
fast_transport::TransportOptions fast_options(bool audit=false) {
 fast_transport::TransportOptions o;
 if(audit){o.guard_extremum_bounds=true;o.geometry_rtol=o.geometry_atol=1e-11;o.quadrature_rtol=1e-9;o.quadrature_atol=1e-11;o.max_step_fraction=.025;o.fast_tail_u=INFINITY;o.fast_split_field_knots=true;}
 return o;
}
transport::Options reference_options(){transport::Options o;o.method=transport::Method::event_mesh;o.resolution=.005;o.tolerance=1e-8;o.probes=16;return o;}
void parity_record(std::ofstream &inputs,std::ofstream &outputs,int id,double b0,Photon const &p,fast_transport::TransportResult const &r,double initial=0,double target=INFINITY) {
 inputs<<id<<' '<<b0;
 for(auto v:{p.n,p.e1,p.e2}) inputs<<' '<<v.x<<' '<<v.y<<' '<<v.z;
 inputs<<' '<<p.r<<' '<<p.psi<<' '<<p.alpha<<' '<<initial<<' '<<p.omega_inf<<' '<<int(p.pol)<<' '<<(std::isinf(target)?-1:target)<<'\n';
 outputs<<id<<','<<r.tau<<','<<r.distance<<','<<r.state[0]<<','<<r.state[1]<<','<<r.state[2]<<','<<int(r.termination)<<','<<r.stats.geometry_steps<<','<<r.stats.field_evaluations<<','<<r.stats.quadrature_intervals<<'\n';
}
}
int main() {
 try {
 BField field("table/bfield_t10.txt",B_pole,R_star);
 std::ofstream inputs("output/fast_parity_inputs.txt"),ported("output/fast_parity_port.csv");
 inputs<<std::setprecision(17);ported<<std::setprecision(17)<<"id,tau,distance,r,psi,alpha,outcome,steps,fields,quadrature_intervals\n";
 int parity_id=0,failures=0;
 std::ifstream table("table/bench_od.txt");std::string line;
 while(std::getline(table,line)) {
  if(line.empty()||line[0]=='#')continue;
  Case c;std::istringstream row(line);if(!(row>>c.b0>>c.muz>>c.pol>>c.energy>>c.reference))throw std::runtime_error("Bad radial reference");
  auto p=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol?Polarization::E:Polarization::O);
  auto result=fast_transport::propagate(field,Boltzmann(c.b0),p);
  parity_record(inputs,ported,parity_id++,c.b0,p,result);
 }
 std::vector<AngularCase> cases;
 std::mt19937_64 rng(20260924);std::uniform_real_distribution<double> U(0,1);
 for(int i=0;i<80;++i) {
  double b0=std::array<double,4>{-.1,-.5,-.9,.3}[i%4];
  double muz=2*U(rng)-1,alpha=pi*U(rng),az=2*pi*U(rng),r=30,energy=std::pow(10.,-1+3*U(rng));
  if(i<16){muz=std::array<double,4>{-1.,-.5,0.,1.}[i%4];alpha=std::array<double,4>{0.,pi/2,2.7,3.05}[i/4];az=.7;energy=1;}
  if(i>=72){r=R_star;alpha=std::array<double,4>{.2,.6,pi/2,3.05}[i%4];}
  for(auto pol:{Polarization::E,Polarization::O})cases.push_back({i,b0,muz,az,transport::make_photon(r,muz,alpha,az,energy,pol)});
 }
 // A restarted photon can have nonzero psi in its retained plane basis.
 for(int i=0;i<6;++i)for(auto pol:{Polarization::E,Polarization::O}) {
  auto p=transport::make_photon(30,.2,.3+.48*i,.7,.3,pol);p.psi=i%2 ? -1.2 : .37;
  cases.push_back({80+i,i%2 ? -.3 : .3,.2,.7,p});
 }
 // Rays just above/below stellar grazing, plus exactly inward radial rays.
 for(int i=0;i<4;++i)for(auto pol:{Polarization::E,Polarization::O}) {
  double peri=R_star*(1+(i%2 ? 1e-6 : -1e-6));
  double impact=peri/std::sqrt(1-rs/peri);
  double alpha=i<2 ? pi-std::asin(impact*std::sqrt(1-rs/30)/30) : pi;
  cases.push_back({86+i,-.1,i==3 ? -1. : .2,.7,transport::make_photon(30,i==3 ? -1. : .2,alpha,.7,1,pol)});
 }
 // Regression for the historical unbracketed angular extremum at the pole.
 int guarded_regressions=0;
 for(auto const &c:cases)if(c.id==89) {
  auto o=fast_options(true);o.guard_extremum_bounds=false;bool observed=false;
  try{fast_transport::propagate(field,Boltzmann(c.b0),c.photon,o);}
  catch(std::exception const &e){observed=std::string(e.what())=="Root must be bracketed";}
  o.guard_extremum_bounds=true;
  auto fixed=fast_transport::propagate(field,Boltzmann(c.b0),c.photon,o);
  if(!observed || fixed.termination!=fast_transport::TransportTermination::absorbed || fixed.tau!=0 || std::abs(fixed.state[0]-R_star)>1e-9)throw std::runtime_error("Angular extremum bounds regression failed");
  ++guarded_regressions;
 }
 std::map<double,Edges> edges;
 for(auto const &c:cases)if(!edges.count(c.b0))edges[c.b0]=distribution_edges(Boltzmann(c.b0));
 std::vector<transport::Result> truth;
 for(auto const &c:cases)truth.push_back(transport::integrate(field,Boltzmann(c.b0),c.photon,reference_options()));
 std::cout<<"Fine reference completed: "<<cases.size()<<" general cases\n"<<std::flush;
 std::ofstream angular("output/fast_angular_comparison.csv"),timing("output/fast_angular_timing.csv");
 angular<<std::setprecision(17)<<"id,b0,muz,alpha,azimuth,r,psi,energy,pol,method,tau,reference,abs_error,relative_error,microseconds,steps,fields,outcome,reference_outcome,status\n";
 timing<<std::setprecision(17)<<"method,cases,batch_median_ms,batch_min_ms,batch_max_ms\n";
 for(std::string method:{"event_guard","phase_cap","ft07","transport_fast","transport_fast_guarded","transport_fast_audit"}) {
  auto raw_solve=[&](AngularCase const &c) {
   Boltzmann fb(c.b0);
   if(method=="transport_fast" || method=="transport_fast_guarded" || method=="transport_fast_audit") {
    auto options=fast_options(method=="transport_fast_audit");
    if(method=="transport_fast_guarded")options.guard_extremum_bounds=true;
    auto r=fast_transport::propagate(field,fb,c.photon,options);
    return Answer{r.tau,r.distance,int(r.termination),long(r.stats.field_evaluations),long(r.stats.geometry_steps)};
   }
   transport::Options o;o.method=method=="event_guard"?transport::Method::event_guard:method=="phase_cap"?transport::Method::phase_cap:transport::Method::ft07;
   o.resolution=method=="ft07"?.01:.1;
   auto r=transport::integrate(field,fb,c.photon,o,&edges.at(c.b0));
   return Answer{r.tau,r.length,int(r.outcome),r.evaluations,r.steps};
  };
  auto solve=[&](AngularCase const &c){try{return raw_solve(c);}catch(std::exception const &e){return Answer{NAN,NAN,-1,0,0,e.what()};}};
  std::vector<Answer> results;std::vector<std::vector<double>> times(cases.size());std::vector<double> batches;
  for(auto const &c:cases)results.push_back(solve(c)); // warmup
  for(int rep=0;rep<5;++rep) {
   auto batch=Clock::now();
   for(size_t i=0;i<cases.size();++i){auto start=Clock::now();auto r=solve(cases[i]);times[i].push_back(std::chrono::duration<double,std::micro>(Clock::now()-start).count());if(r.status!=results[i].status || (r.status=="ok" && (r.tau!=results[i].tau||r.outcome!=results[i].outcome)))throw std::runtime_error("Non-deterministic angular solve");}
   batches.push_back(std::chrono::duration<double,std::milli>(Clock::now()-batch).count());
  }
  std::sort(batches.begin(),batches.end());timing<<method<<','<<cases.size()<<','<<batches[2]<<','<<batches.front()<<','<<batches.back()<<'\n';timing.flush();
  for(size_t i=0;i<cases.size();++i) {
   auto const &c=cases[i];auto const &p=c.photon;auto const &r=results[i];std::sort(times[i].begin(),times[i].end());
   double ae=std::abs(r.tau-truth[i].tau),re=ae/std::max(1e-12,std::abs(truth[i].tau));
   // FT07 truncates particle tails; report its full-domain bias without hiding it.
   bool ok=r.status=="ok" && r.outcome==int(truth[i].outcome)&&(method=="ft07" || ae<=1e-7+1e-4*std::abs(truth[i].tau));
   failures+=!ok;
   angular<<c.id<<','<<c.b0<<','<<c.muz<<','<<p.alpha<<','<<c.azimuth<<','<<p.r<<','<<p.psi<<','<<p.omega_inf<<','<<int(p.pol)<<','<<method<<','<<r.tau<<','<<truth[i].tau<<','<<ae<<','<<re<<','<<times[i][2]<<','<<r.steps<<','<<r.fields<<','<<r.outcome<<','<<int(truth[i].outcome)<<','<<(r.status!="ok"?r.status:ok?"ok":"mismatch")<<'\n';
   if(method=="transport_fast")parity_record(inputs,ported,parity_id++,c.b0,p,fast_transport::propagate(field,Boltzmann(c.b0),p));
   if(!ok)std::cerr<<"Angular mismatch "<<c.id<<' '<<int(p.pol)<<' '<<method<<" abs="<<ae<<" rel="<<re<<" status="<<r.status<<" outcome="<<r.outcome<<'/'<<int(truth[i].outcome)<<'\n';
  }
  angular.flush();std::cout<<method<<" general batch_ms="<<batches[2]<<" failures="<<failures<<'\n'<<std::flush;
 }
 std::ofstream scatter("output/fast_scattering_validation.csv");
 scatter<<std::setprecision(17)<<"id,U,initial_tau,target,total_tau,scatters,distance,r,reference_cumulative,absolute_residual,reported_tau_error,impact_relative_error,status\n";
 int scattering=0,escapes=0;
 for(int i=0;i<12;++i) {
  double b0=i%3==0 ? .3 : -.3;Boltzmann fb(b0);
  auto p=transport::make_photon(30,i%2 ? -.4 : .4,.25+.24*i,.7+.31*i,.3,i%2 ? Polarization::E : Polarization::O);
  auto ref=reference_options();auto total=transport::integrate(field,fb,p,ref);
  double initial=.37;
  for(double u:{.1,.5,.9}) {
   double target=initial-std::log(u);
   auto r=fast_transport::propagate(field,fb,p,{},target,initial);
   parity_record(inputs,ported,parity_id++,b0,p,r,initial,target);
   bool expected=target<initial+total.tau,actual=r.termination==fast_transport::TransportTermination::scattered;
   double cumulative=total.tau;
   if(actual){auto limited=ref;limited.path_limit=r.distance;cumulative=transport::integrate(field,fb,p,limited).tau;++scattering;}else ++escapes;
   double residual=actual ? std::abs(initial+cumulative-target) : std::abs(initial+total.tau-r.tau);
   double b1=p.r*std::sin(p.alpha)/std::sqrt(1-rs/p.r),b2=r.state[0]*std::sin(r.state[2])/std::sqrt(1-rs/r.state[0]);
   double invariant=std::abs(b1-b2)/b1;
   bool ok=expected==actual && residual<1e-6+1e-5*target && invariant<1e-7;
   failures+=!ok;
   scatter<<i<<','<<u<<','<<initial<<','<<target<<','<<total.tau<<','<<actual<<','<<r.distance<<','<<r.state[0]<<','<<cumulative<<','<<residual<<','<<(actual?std::abs(r.tau-target):0)<<','<<invariant<<','<<(ok?"ok":"mismatch")<<'\n';
  }
 }
 // An already reached absolute target must return the starting state unchanged.
 auto p=cases[50].photon;auto r=fast_transport::propagate(field,Boltzmann(cases[50].b0),p,{},.2,.37);
 parity_record(inputs,ported,parity_id++,cases[50].b0,p,r,.37,.2);
 if(r.distance!=0 || r.tau!=.37 || r.termination!=fast_transport::TransportTermination::scattered)++failures;
 std::ofstream summary("output/fast_validation_summary.json");summary<<"{\"general_cases\": "<<cases.size()<<", \"parity_cases\": "<<parity_id<<", \"scatter_events\": "<<scattering<<", \"no_scatter_events\": "<<escapes<<", \"guarded_regressions\": "<<guarded_regressions<<", \"failures\": "<<failures<<"}\n";
 std::cout<<"Validation completed: failures="<<failures<<", parity="<<parity_id<<'\n';
 return failures?1:0;
 }catch(std::exception const &e){std::cerr<<e.what()<<'\n';return 1;}
}
