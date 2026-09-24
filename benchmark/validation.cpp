#include "../output/velocity_integral.hpp"
#include "../src/radial_optical_depth.hpp"
#include <cassert>

void require(bool condition,char const *message) { if(!condition) throw std::runtime_error(message); }

int main() {
 try {
  BField field("table/bfield_t10.txt",B_pole,R_star);
  std::ofstream ft("output/ft07_truncation.csv"), cumulative("output/cumulative.csv"), events("output/events.csv"), stress("output/stress.csv");
  for(auto *out:{&ft,&cumulative,&events,&stress}) { require(bool(*out),"Cannot open validation output"); *out<<std::setprecision(17); }
  ft<<"b0,muz,pol,omega_inf,full,truncated,omitted_fraction,low,median,high\n";
  cumulative<<"b0,muz,pol,omega_inf,r,spatial,velocity,rel_error\n";
  events<<"b0,muz,pol,omega_inf,U,tau_draw,scatters,r_spatial,r_velocity,relative_radius_error\n";
  stress<<"b0,muz,pol,omega_inf,spatial,velocity,rel_error\n";
  std::ifstream in("table/bench_od.txt"); std::string line;
  std::map<double,Edges> edges;
  double max_trunc=0,max_cumulative=0,max_event=0,max_stress=0,max_internal=0;
  int n_full=0,n_cumulative=0,n_events=0,n_stress=0,n_escape=0;
  while(std::getline(in,line)) {
    if(line.empty()||line[0]=='#')continue;
    Case c;std::istringstream row(line); row>>c.b0>>c.muz>>c.pol>>c.energy>>c.reference;
    Boltzmann fb(c.b0);
    if(!edges.count(c.b0))edges[c.b0]=distribution_edges(fb);
    auto e=edges.at(c.b0);
    RadialOpticalDepth s(field,fb,c.muz,c.energy,c.pol ? Polarization::E : Polarization::O);
    double full=velocity_integral(s.ray,1e-11).tau;
    double spatial=s.integrate_to(10000,1e-12,1e-11);
    double internal=std::abs(spatial-full)/full;
    max_internal=std::max(max_internal,internal);
    require(internal<1e-8,"Total-depth cross-check failed");
    double truncated=velocity_integral(s.ray,1e-11,10000,e.low,e.high).tau;
    double omitted=(full-truncated)/full;
    max_trunc=std::max(max_trunc,omitted);++n_full;
    ft<<c.b0<<','<<c.muz<<','<<c.pol<<','<<c.energy<<','<<full<<','<<truncated<<','<<omitted<<','<<e.low<<','<<e.median<<','<<e.high<<'\n';
    if(!(c.b0==-.1 || c.b0==-.5 || c.b0==-.9) || !(c.muz==0 || c.muz==.4 || c.muz==.9) || c.energy==.1 || c.energy==10)continue;
    for(double k:{3.,1.,.3,.03}) {
      double beta=std::max(s.ray.beta(R_star)*.999,k*c.b0);
      double r=s.ray.radius_for_x(s.ray.g(beta));
      double a=s.integrate_to(r,1e-12,1e-10),b=velocity_integral(s.ray,1e-10,r).tau;
      double error=std::abs(a-b)/std::max(1e-12,b);
      max_cumulative=std::max(max_cumulative,error);++n_cumulative;
      require(std::abs(a-b)<1e-9+1e-8*b,"Cumulative-depth cross-check failed");
      cumulative<<c.b0<<','<<c.muz<<','<<c.pol<<','<<c.energy<<','<<r<<','<<a<<','<<b<<','<<error<<'\n';
    }
    for(double U:{.01,.1,.5,.9}) {
      double draw=-std::log(U),r1=0,r2=0,error=0;
      bool scatters=draw<full;
      require(scatters==(draw<spatial),"Escape decision differs");
      if(scatters) {
        r1=zriddr([&](double r){return s.integrate_to(r,1e-11,1e-10)-draw;},R_star,s.r_end,1e-9);
        r2=zriddr([&](double r){return velocity_integral(s.ray,1e-10,r).tau-draw;},R_star,s.r_end,1e-9);
        error=std::abs(r1-r2)/r2;max_event=std::max(max_event,error);
        require(error<1e-8,"Event-radius cross-check failed");
      } else ++n_escape;
      ++n_events;
      events<<c.b0<<','<<c.muz<<','<<c.pol<<','<<c.energy<<','<<U<<','<<draw<<','<<scatters<<','<<r1<<','<<r2<<','<<error<<'\n';
    }
  }
  require(n_full==900,"Validation did not read all reference cases");
  for(double b0:{-.03,-.05,-.99})for(double muz:{0.,1e-6,.01,.9})for(double energy:{.001,1.,1000.})for(auto pol:{Polarization::E,Polarization::O}) {
    Boltzmann fb(b0);RadialOpticalDepth s(field,fb,muz,energy,pol);
    double a=s.integrate_to(10000,1e-12,1e-10),b=velocity_integral(s.ray,1e-10).tau;
    double error=std::abs(a-b)/std::max(1e-12,b);max_stress=std::max(max_stress,error);++n_stress;
    require(std::abs(a-b)<1e-9+1e-8*b,"Stress cross-check failed");
    stress<<b0<<','<<muz<<','<<int(pol==Polarization::E)<<','<<energy<<','<<a<<','<<b<<','<<error<<'\n';
  }
  // A zero RHS must still obey the cap, including after error=0 grows h_new.
  auto zero=[](double,std::array<double,1> const &,std::array<double,1> &dy){dy[0]=0;};
  StepperDopr5<1,decltype(zero)> stepper(zero,1e-8,1e-8);
  for(double direction:{1.,-1.}) {
    stepper.init(0,direction*10,{0});
    for(int i=0;i<3;++i) {stepper.do_step(.025); require(std::abs(stepper.h_old)<=.025 && stepper.h_old*direction>0,"Cap not enforced");stepper.update_old();}
  }
  bool caught=false;try{stepper.do_step(0);}catch(std::invalid_argument const &){caught=true;}
  require(caught,"Invalid cap accepted");
  std::ofstream summary("output/validation_summary.txt");
  summary<<std::setprecision(17)<<"all_passed=true\nfull_cases="<<n_full<<"\nmax_spatial_velocity_relative="<<max_internal<<"\nmax_FT07_truncated_fraction="<<max_trunc<<"\ncumulative_cases="<<n_cumulative<<"\nmax_cumulative_relative="<<max_cumulative<<"\nevent_cases="<<n_events<<"\nescapes="<<n_escape<<"\nmax_event_radius_relative="<<max_event<<"\nstress_cases="<<n_stress<<"\nmax_stress_relative="<<max_stress<<"\nstep_cap_tests=passed\n";
  std::cout<<"Validation passed: "<<n_full<<" totals, "<<n_cumulative<<" cumulative, "<<n_events<<" events, "<<n_stress<<" stress cases.\n";
 }catch(std::exception const &e) {std::cerr<<"validation: "<<e.what()<<'\n';return 1;}
}
