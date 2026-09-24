#pragma once
#include "resonance_panels.hpp"
#include "dopr5.hpp"

namespace transport {
enum class Method { ft07, phase_cap, event_mesh, event_guard };
enum class Outcome { escaped, surface, scattered, path_limit };
struct Options {
  Method method=Method::event_guard;
  double tolerance=1e-6;
  double resolution=.1;
  double initial_step=.01;
  double escape_radius=10000;
  double path_limit=std::numeric_limits<double>::infinity();
  double tau_target=std::numeric_limits<double>::infinity();
  int probes=8;
  bool truncate_distribution=false; // only for matching FT07's retained 99.8%
};
struct Result {
  double tau=0,length=0;
  State state{};
  long evaluations=0,steps=0,quadrature_evaluations=0;
  int caustics=0,double_branches=0,turns=0;
  Outcome outcome=Outcome::escaped;
};

inline std::array<double,3> distribution_edges(Boltzmann const &fb) {
  Boltzmann positive(std::abs(fb.b0));
  auto quantile=[&](double probability) {
    return zriddr([&](double b){
      long n=0;
      return quad::adaptive([&](double v){return positive.f(v);},0.,b,1e-13,1e-12,n)-probability;
    },1e-16,1.,1e-14);
  };
  double a=quantile(.001),b=quantile(.5),c=quantile(.999);
  return fb.b0<0 ? std::array<double,3>{-c,-b,-a} : std::array<double,3>{a,b,c};
}

inline double ft07_cap(Medium const &medium,State const &y,std::array<double,3> const &edges,double resolution) {
  auto p=medium.at(y);double cap=y[0]/10,D=p.discriminant();
  if(D<=1e-6)return cap; // caustic neighborhood is handled by event quadrature
  auto rate=medium.rates(y);double dx=rate[0],dm=rate[1];
  for(double b:medium.branches(p,D)) {
    double db=((1-b*p.mu)*dx+b*dm)*(1-b*b)/(b-p.mu);
    if(!std::isfinite(db) || db==0)continue;
    if(b>=edges[0] && b<=edges[2]) {
      double delta=resolution*(1-b*b)*std::abs(b);
      if(delta==0)delta=resolution*std::min(edges[1]-edges[0],edges[2]-edges[1]);
      cap=std::min(cap,delta/std::abs(db));
      if(dm!=db)cap=std::min(cap,resolution*std::abs(p.mu-b)/std::abs(dm-db));
    } else {
      double edge=db>0 ? edges[0] : edges[2];
      if((edge-b)/db>0 && (edge-b)/db<cap) {
        double width=db>0 ? edges[1]-edges[0] : edges[2]-edges[1];
        cap=std::min(cap,(std::abs(edge-b)+resolution*width)/std::abs(db));
      }
    }
  }
  return cap;
}

inline Result integrate(BField const &field,Boltzmann const &fb,Photon const &photon,
                        Options const &options={},std::array<double,3> const *precomputed_edges=nullptr) {
  if(!(options.tolerance>0 && options.resolution>0 && options.initial_step>0 && options.probes>=4))
    throw std::invalid_argument("Invalid transport options");
  if(!(photon.r>=R_star && photon.alpha>=0 && photon.alpha<=pi && photon.omega_inf>0))
    throw std::invalid_argument("Invalid photon state");
  Result result;result.state={photon.r,photon.psi,photon.alpha};
  if(photon.r==R_star && photon.alpha>pi/2){result.outcome=Outcome::surface;return result;}
  if(photon.r>=options.escape_radius && photon.alpha<=pi/2)return result;
  if(options.tau_target==0){result.outcome=Outcome::scattered;return result;}
  Medium medium{field,fb,photon};
  std::array<double,3> edges{};
  bool truncated=options.method==Method::ft07 || options.truncate_distribution;
  if(truncated)edges=precomputed_edges ? *precomputed_edges : distribution_edges(fb);
  double low=truncated ? edges[0] : fb.b_min;
  double high=truncated ? edges[2] : fb.b_max;
  std::vector<double> velocities={0};
  for(double k:{.0625,.125,.25,.5,1.,2.,4.,8.})if(std::abs(k*fb.b0)<1)velocities.push_back(k*fb.b0);
  for(double v:{.9,.99,.999,.9999})velocities.push_back(std::copysign(v,fb.b0));
  if(truncated){velocities.push_back(edges[0]);velocities.push_back(edges[2]);}
  auto deriv=[](double,State const &y,State &dy){dy=geodesic(y);};
  double geometry_tolerance=std::min(1e-10,options.tolerance*.01);
  StepperDopr5<3,decltype(deriv)> stepper(deriv,geometry_tolerance,geometry_tolerance);
  stepper.add_event([&](double,State const &y){return y[0]-options.escape_radius;});
  stepper.add_event([](double,State const &y){return y[0]-R_star;});
  stepper.init(0,options.initial_step,result.state);
  while(true) {
    if(++result.steps>2000000)throw std::runtime_error("Transport step budget exceeded");
    double r=stepper.y_old[0],cap=options.resolution*r;
    if(options.method==Method::event_guard)cap=std::min(.4,options.resolution*4)*r;
    if(options.method==Method::ft07)cap=ft07_cap(medium,stepper.y_old,edges,options.resolution);
    if(options.method==Method::phase_cap) {
      auto rates=medium.rates(stepper.y_old);
      // For all angles, resolve changes of BOTH x and mu on the thermal
      // scale, with a geometry cap protecting stationary derivatives.
      double speed=std::abs(rates[0])+std::abs(rates[1])/std::max(.01,std::abs(fb.b0));
      auto p=medium.at(stepper.y_old);
      double fmin=1e100,fmax=-1e100;
      for(double beta:velocities) {
        if(std::abs(beta)>std::min(.999,4*std::abs(fb.b0)))continue;
        double f=p.surface(beta);fmin=std::min(fmin,f);fmax=std::max(fmax,f);
      }
      double distance=fmin>0 ? fmin : fmax<0 ? -fmax : 0;
      cap=std::min(r/10,(options.resolution*fb.b0*fb.b0+.5*distance)/std::max(speed,1e-30));
    }
    cap=std::min({cap,.4*r,options.path_limit-stepper.x_old});
    if(cap<=0){result.outcome=Outcome::path_limit;break;}
    stepper.do_step(cap);
    int event=stepper.detect_event();
    // detect_event may shorten h_old while preserving the original dense
    // polynomial. Do not prepare it again after event truncation.
    if(event<0)stepper.prepare_dense();
    // An inward grazing ray can hit the star and turn within one accepted
    // step even if BOTH radii exceed R_star. Check the radial minimum too.
    if(stepper.y_old[2]>pi/2 && stepper.y_new[2]<=pi/2) {
      double turn=zriddr([&](double l){return stepper.dense_out(l)[2]-pi/2;},stepper.x_old,stepper.x_old+stepper.h_old,1e-12);
      if(stepper.dense_out(turn)[0]<=R_star) {
        double hit=zriddr([&](double l){return stepper.dense_out(l)[0]-R_star;},stepper.x_old,turn,1e-12);
        stepper.h_old=hit-stepper.x_old;stepper.y_new=stepper.dense_out(hit);event=1;
      }
    }
    // Crossing the outer boundary inward is entry, not escape. The dense
    // polynomial already represents the shortened step and is kept intact.
    if(event==0 && stepper.y_new[2]>pi/2)event=-2;
    double h=stepper.h_old;
    auto path=[&](double t){return stepper.dense_out(stepper.x_old+t*h);};
    ResonancePanels panels(medium,path,h,velocities,options.probes,low,high);
    for(size_t i=1;i<panels.boundaries.size();++i) {
      auto a=panels.boundaries[i-1],b=panels.boundaries[i];
      // Allocate absolute error by proper-path fraction at the local radius.
      double atol=options.tolerance*.01*std::min(1.,h/r)/(panels.boundaries.size()-1);
      double dtau=panels.integrate(a,b,atol,options.tolerance);
      if(result.tau+dtau>=options.tau_target) {
        double needed=options.tau_target-result.tau;
        double t=zriddr([&](double t){return t==b.t ? dtau-needed : panels.integrate(a,{t,false},atol,options.tolerance)-needed;},a.t,b.t,1e-11);
        result.tau=options.tau_target;result.length=stepper.x_old+t*h;result.state=path(t);
        result.outcome=Outcome::scattered;result.evaluations=medium.evaluations;
        result.quadrature_evaluations+=panels.evaluations;return result;
      }
      result.tau+=dtau;
    }
    result.caustics+=panels.caustics;result.double_branches+=panels.double_branches;
    result.quadrature_evaluations+=panels.evaluations;
    result.turns+=(stepper.y_old[2]>pi/2 && stepper.y_new[2]<=pi/2);
    result.length=stepper.x_old+h;result.state=stepper.y_new;
    if(event>=0){result.outcome=event==0 ? Outcome::escaped : Outcome::surface;break;}
    if(result.length>=options.path_limit){result.outcome=Outcome::path_limit;break;}
    stepper.update_old();
  }
  result.evaluations=medium.evaluations;
  return result;
}
}
