#include "../output/benchmark_common.hpp"
#include <random>

void check(bool value,char const *message){if(!value)throw std::runtime_error(message);}

// Synthetic medium with a known narrow allowed pocket between a and b.
// Geometry and all external probes may be outside that pocket. The exact
// support is used ONLY by the independent analytic-variable reference.
struct Pocket : transport::Medium {
  double centre,width;
  transport::ResonancePoint at(transport::State const &y) const {
    double d=y[0]-centre,D=.5*(.25*width*width-d*d);
    return {.5*std::log1p(D),0,1};
  }
};

int main() {
 try {
  BField field("table/bfield_t10.txt",B_pole,R_star);
  std::ofstream pocket_file("output/pocket_validation.csv"),event_file("output/scattering_validation.csv");
  pocket_file<<std::setprecision(17)<<"width,numerical,independent_reference,relative_error,detected_caustics\n";
  event_file<<std::setprecision(17)<<"id,U,scatters,tau_total,tau_target,length,r,reference_cumulative,error\n";
  Boltzmann fb(-.1);auto photon=transport::make_photon(30,.4,.6,.8,1,Polarization::E);
  double max_pocket=0,max_opacity=0,max_event=0,max_invariant=0;
  for(double width:{.1,.001,.00001,.000001}) {
    Pocket material{{field,fb,photon},.371239,width};
    auto path=[](double t){return transport::State{t,0,0};};
    transport::ResonancePanels panels(material,path,1.,std::vector<double>{0,-.05,-.1,-.2},8);
    double numerical=0;
    for(size_t i=1;i<panels.boundaries.size();++i)numerical+=panels.integrate(panels.boundaries[i-1],panels.boundaries[i],1e-9,1e-7);
    long calls=0;
    double reference=quadrature::adaptive([&](double theta){
      double D=.125*width*width*std::cos(theta)*std::cos(theta);
      transport::ResonancePoint p{.5*std::log1p(D),0,1};
      return material.opacity(p,D)*.5*width*std::cos(theta);
    },-pi/2,pi/2,1e-11,1e-10,calls);
    double error=std::abs(numerical-reference)/reference;max_pocket=std::max(max_pocket,error);
    pocket_file<<width<<','<<numerical<<','<<reference<<','<<error<<','<<panels.caustics<<'\n';pocket_file.flush();
    check(panels.caustics==2,"Failed to discover same-sign pocket boundaries");
    check(error<3e-5,"Narrow-pocket optical depth mismatch");
  }
  // Compare the new stable branch sum to the original note/source formula
  // at arbitrary angles away from an exact multiple root.
  std::mt19937_64 rng(1707);std::uniform_real_distribution<double> U(0,1);
  for(int i=0;i<3000;++i) {
    double b0=i%2 ? -.5 : .5;Boltzmann fb(b0);
    auto p=transport::make_photon(10+200*U(rng),2*U(rng)-1,pi*U(rng),2*pi*U(rng),.01*std::pow(10.,4*U(rng)),i%3 ? Polarization::E : Polarization::O);
    transport::Medium medium{field,fb,p};auto rp=medium.at({p.r,p.psi,p.alpha});
    if(std::abs(rp.discriminant())<1e-7)continue;
    PhotonEvolution original{field,fb,p.n,p.e1,p.e2,{p.r,p.psi,p.alpha,0},p.omega_inf,p.pol};
    double a=medium.opacity(rp,rp.discriminant()),b=original.calc_dtaudl(p.n,p.e1,p.e2,p.r,p.psi,p.alpha,p.omega_inf,p.pol);
    double error=std::abs(a-b)/std::max(1e-12,std::abs(b));max_opacity=std::max(max_opacity,error);
    check(std::abs(a-b)<1e-11+1e-8*std::abs(b),"General-angle opacity formula changed");
  }
  int event_count=0,escapes=0;
  for(int i=0;i<12;++i) {
    Boltzmann population(i%3==0 ? .3 : -.3);
    auto p=transport::make_photon(30,i%2 ? -.4 : .4,.25+.24*i,.7+.31*i,.3,i%2 ? Polarization::E : Polarization::O);
    transport::Options options;options.tolerance=1e-7;
    auto total=transport::integrate(field,population,p,options);
    double impact_initial=p.r*std::sin(p.alpha)/std::sqrt(1-rs/p.r);
    double impact_final=total.state[0]*std::sin(total.state[2])/std::sqrt(1-rs/total.state[0]);
    double invariant=std::abs(impact_initial-impact_final)/impact_initial;max_invariant=std::max(max_invariant,invariant);
    check(invariant<1e-7,"Geodesic impact parameter not conserved");
    for(double u:{.1,.5,.9}) {
      ++event_count;double target=-std::log(u);
      auto cfg=options;cfg.tau_target=target;
      auto event=transport::integrate(field,population,p,cfg);
      bool scatters=target<total.tau;
      check(scatters==(event.outcome==transport::Outcome::scattered),"Incorrect scatter/escape decision");
      double ref=0,error=0;
      if(scatters) {
        cfg=options;cfg.method=transport::Method::event_mesh;cfg.resolution=.005;cfg.probes=16;cfg.tolerance=1e-8;cfg.path_limit=event.length;
        ref=transport::integrate(field,population,p,cfg).tau;
        error=std::abs(ref-target);max_event=std::max(max_event,error);
        check(error<1e-6+1e-5*target,"Scattering event cumulative-depth mismatch");
      }else ++escapes;
      event_file<<i<<','<<u<<','<<scatters<<','<<total.tau<<','<<target<<','<<event.length<<','<<event.state[0]<<','<<ref<<','<<error<<'\n';event_file.flush();
    }
  }
  // Surface interception must detect an inward grazing ray even if a
  // proposed step straddles periapsis with both endpoints above the star.
  for(double offset:{-1e-6,1e-6}) {
    double peri=R_star*(1+offset),impact=peri/std::sqrt(1-rs/peri);
    double alpha=pi-std::asin(impact*std::sqrt(1-rs/30)/30);
    auto grazing=transport::make_photon(30,.2,alpha,.7,1,Polarization::E);
    auto answer=transport::integrate(field,fb,grazing);
    check(answer.outcome==(offset<0 ? transport::Outcome::surface : transport::Outcome::escaped),"Grazing surface interception failed");
  }
  // A full-domain comparison isolates FT07's particle-tail truncation using
  // the SAME general-angle solver with a different controller, never the
  // deleted radial velocity inversion.
  std::ofstream decomp("output/ft07_decomposition.csv");decomp<<std::setprecision(17);
  decomp<<"b0,muz,pol,omega_inf,table_reference,full,truncated,omitted_fraction\n";
  std::ifstream input("table/bench_od.txt");std::string line;int decomp_count=0;
  std::map<double,std::array<double,3>> edge_cache;
  while(std::getline(input,line)) {
    if(line.empty()||line[0]=='#')continue;
    Case c;std::istringstream row(line);row>>c.b0>>c.muz>>c.pol>>c.energy>>c.reference;
    Boltzmann population(c.b0);if(!edge_cache.count(c.b0))edge_cache[c.b0]=transport::distribution_edges(population);
    auto p=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
    transport::Options cfg;cfg.tolerance=1e-8;
    double full=transport::integrate(field,population,p,cfg).tau;
    cfg.truncate_distribution=true;
    double truncated=transport::integrate(field,population,p,cfg,&edge_cache.at(c.b0)).tau;
    decomp<<c.b0<<','<<c.muz<<','<<c.pol<<','<<c.energy<<','<<c.reference<<','<<full<<','<<truncated<<','<<(full-truncated)/full<<'\n';++decomp_count;
  }
  check(decomp_count==900,"Incomplete FT07 decomposition");
  std::ofstream summary("output/transport_validation_summary.json");
  summary<<std::setprecision(17)<<"{\n\"passed\": true,\n\"max_pocket_relative_error\": "<<max_pocket<<",\n\"max_opacity_relative_difference\": "<<max_opacity<<",\n\"max_event_absolute_error\": "<<max_event<<",\n\"max_impact_invariant_relative_error\": "<<max_invariant<<",\n\"event_tests\": "<<event_count<<",\n\"escape_tests\": "<<escapes<<",\n\"grazing_surface_tests\": 2,\n\"ft07_decomposition_cases\": "<<decomp_count<<"\n}\n";
  std::cout<<"Transport validation passed\n";
 }catch(std::exception const &e){std::cerr<<"transport validation: "<<e.what()<<'\n';return 1;}
}
