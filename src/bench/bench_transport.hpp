#pragma once
#include "bench_panels.hpp"
#include "../dopr5.hpp"

namespace transport {
enum class Method { ft07, geometric, panels, reference, baseline };
struct Options {
  double tolerance=1e-6, resolution=.2, initial_step=.01;
  double escape_radius=10000, path_limit=INFINITY, tau_target=INFINITY;
  int probes=8;
};
struct Result {
  double tau=0,length=0;
  State state{};
  long evaluations=0,opacity_evaluations=0,steps=0,quadrature_evaluations=0;
  long caustics=0,double_branches=0,turns=0;
  int outcome=0; // 0 escape, 1 star, 2 path limit, 3 scattering
};

// FT07 eqs. (31), (38), (39), with d ln x/dl replacing d ln B/dl
// to include gravitational redshift. dmu/dl includes bending and field rotation.
// The fiducial 1/100 is multiplied by resolution/.01 for convergence tests.
// Distribution bands are supplied externally and may describe several peaks.
inline double ft07_cap(Medium const &medium,State const &y,double resolution) {
  auto p=medium.at(y); double cap=y[0]/10, D=p.discriminant();
  // A simple caustic is reached only after infinitely many relative-gap
  // steps. Within this small neighborhood hand it to the event quadrature
  // (sin^2 map), rather than discarding the finite endpoint optical depth.
  if(D<=1e-8)return cap;
  auto [dx,dm]=medium.rates(y);
  for(double beta:medium.branches(p,D)) {
    double db=((1-beta*p.mu)*dx+beta*dm)*(1-beta*beta)/(beta-p.mu);
    if(!std::isfinite(db) || db==0)continue;
    bool active=false;
    for(auto band:medium.distribution.bands) {
      if(beta>=band[0] && beta<=band[2]) {
        active=true;
        double delta=resolution*(1-beta*beta)*std::abs(beta);
        if(delta==0)delta=resolution*std::min(band[1]-band[0],band[2]-band[1]);
        cap=std::min(cap,delta/std::abs(db));
      } else {
        double edge=db>0 ? band[0] : band[2];
        double distance=(edge-beta)/db;
        if(distance>0 && distance<cap) {
          double width=db>0 ? band[1]-band[0] : band[2]-band[1];
          cap=std::min(cap,(std::abs(edge-beta)+resolution*width)/std::abs(db));
        }
      }
    }
    if(active && dm!=db)cap=std::min(cap,resolution*std::abs(p.mu-beta)/std::abs(dm-db));
  }
  return cap;
}

// Direct QAGS over each smooth trajectory interval: the algorithm in ref.cpp.
// This wrapper reports errors to the benchmark, instead of aborting the suite.
template<class F> double checked_qags(F &&f,double a,double b,double atol,double rtol,long &evals) {
  auto fn=[&](double x){++evals;return f(x);};
  gsl_function g;
  g.function=+[](double x,void *p){return (*static_cast<decltype(fn)*>(p))(x);};g.params=&fn;
  auto *w=gsl_integration_workspace_alloc(2000);
  if(!w)throw std::runtime_error("GSL workspace allocation failed");
  double result,error;
  int status=gsl_integration_qags(&g,a,b,atol,rtol,2000,w,&result,&error);
  gsl_integration_workspace_free(w);
  if(status || !std::isfinite(result))
    throw std::runtime_error(std::format("QAGS: {} result={:.6e} error={:.6e}",gsl_strerror(status),result,error));
  return result;
}

inline Result integrate(BField const &field,Distribution const &distribution,Photon const &photon,
                        Method method,Options const &o={}) {
  if(!(o.tolerance>0 && o.resolution>0 && o.initial_step>0 && o.probes>=4))
    throw std::invalid_argument("Invalid integration options");
  if(!(std::isfinite(distribution.mean) && distribution.mean!=0))
    throw std::invalid_argument("Current-normalized opacity requires a nonzero mean velocity");
  if(!(photon.r>=R_star && photon.alpha>=0 && photon.alpha<=pi && photon.omega_inf>0 &&
       o.path_limit>0 && o.tau_target>0 && o.escape_radius>photon.r))
    throw std::invalid_argument("Invalid photon or termination bounds");
  Result result;result.state={photon.r,photon.psi,photon.alpha};
  if(photon.r==R_star && photon.alpha>pi/2){result.outcome=1;return result;}
  Medium medium{field,distribution,photon};
  auto finish=[&]() {
    result.evaluations=medium.evaluations;
    result.opacity_evaluations=medium.opacity_evaluations;
    return result;
  };
  if(method==Method::geometric || method==Method::baseline) {
    using State4=std::array<double,4>;
    auto rhs=[&](double,State4 const &y,State4 &dy){
      State z={y[0],y[1],y[2]};auto g=geodesic(z);auto p=medium.at(z);
      dy={g[0],g[1],g[2],medium.opacity(p,p.discriminant())};
    };
    StepperDopr5<4,decltype(rhs)> s(rhs,o.tolerance,o.tolerance);
    s.add_event([&](double,State4 const &y){return y[0]-o.escape_radius;});
    s.add_event([](double,State4 const &y){return y[0]-R_star;});
    if(std::isfinite(o.tau_target))s.add_event([&](double,State4 const &y){return y[3]-o.tau_target;});
    s.init(0,o.initial_step,{photon.r,photon.psi,photon.alpha,0});
    while(true) {
      if(++result.steps>2000000)throw std::runtime_error("RK step budget exceeded");
      double cap=method==Method::geometric ? o.resolution*s.y_old[0] : INFINITY;
      s.do_step(std::min(cap,o.path_limit-s.x_old));
      int event=s.detect_event();
      result.tau=s.y_new[3];result.length=s.x_old+s.h_old;
      result.state={s.y_new[0],s.y_new[1],s.y_new[2]};
      result.turns+=(s.y_old[2]>pi/2 && s.y_new[2]<=pi/2);
      if(event>=0){result.outcome=event==2 ? 3 : event;break;}
      if(result.length>=o.path_limit){result.outcome=2;break;}
      s.update_old();
    }
    return finish();
  }
  std::vector<double> velocities=distribution.knots;
  // The PDF adapter may supply a separate partition for the raw QAGS reference.
  if(method==Method::reference && !distribution.reference_knots.empty())
    velocities=distribution.reference_knots;
  if(distribution.low>-1)velocities.push_back(distribution.low);
  if(distribution.high<1)velocities.push_back(distribution.high);
  velocities.push_back(0);
  for(double v:velocities)if(!std::isfinite(v) || v < -1 || v > 1)
    throw std::invalid_argument("Velocity landmarks must lie in [-1,1]");
  std::erase_if(velocities,[](double v){return std::abs(v)==1;});
  std::sort(velocities.begin(),velocities.end());
  velocities.erase(std::unique(velocities.begin(),velocities.end()),velocities.end());
  auto rhs=[](double,State const &y,State &dy){dy=geodesic(y);};
  double orbit_tol=std::min(1e-10,o.tolerance*.01);
  StepperDopr5<3,decltype(rhs)> s(rhs,orbit_tol,orbit_tol);
  s.add_event([&](double,State const &y){return y[0]-o.escape_radius;});
  s.add_event([](double,State const &y){return y[0]-R_star;});
  s.init(0,o.initial_step,result.state);
  while(true) {
    if(++result.steps>2000000)throw std::runtime_error("Orbit step budget exceeded");
    double r=s.y_old[0];
    double cap=method==Method::ft07 ? ft07_cap(medium,s.y_old,o.resolution) : o.resolution*r;
    cap=std::min({cap,.4*r,o.path_limit-s.x_old});
    if(cap<=0){result.outcome=2;break;}
    s.do_step(cap);
    int event=s.detect_event();
    if(event<0)s.prepare_dense(); // detect_event already preserves the full polynomial
    // Detect an occultation even if the radial minimum is inside a step.
    if(s.y_old[2]>pi/2 && s.y_new[2]<=pi/2) {
      double turn=zriddr([&](double l){return s.dense_out(l)[2]-pi/2;},s.x_old,s.x_old+s.h_old,1e-12);
      if(s.dense_out(turn)[0]<=R_star) {
        double hit=zriddr([&](double l){return s.dense_out(l)[0]-R_star;},s.x_old,turn,1e-12);
        s.h_old=hit-s.x_old;s.y_new=s.dense_out(hit);event=1;
      }
    }
    double h=s.h_old;
    auto path=[&](double t){return s.dense_out(s.x_old+t*h);};
    ResonancePanels panels(medium,path,h,velocities,o.probes,distribution.low,distribution.high);
    for(size_t i=1;i<panels.boundaries.size();++i) {
      auto left=panels.boundaries[i-1],right=panels.boundaries[i];
      double atol=o.tolerance*.01*std::min(1.,h/r)/(panels.boundaries.size()-1);
      auto integrate_panel=[&](Boundary end) {
        if(method!=Method::reference)return panels.integrate(left,end,atol,o.tolerance);
        // Same spatial QAGS integrand as ref.cpp; generalize its radial-only
        // beta=0 event and do not terminate at D=0 (a ray may reenter).
        return checked_qags([&](double t){auto p=medium.at(path(t));return medium.opacity(p,p.discriminant())*h;},
                            left.t,end.t,o.tolerance,o.tolerance,result.quadrature_evaluations);
      };
      double dtau=integrate_panel(right);
      if(result.tau+dtau>=o.tau_target) {
        double needed=o.tau_target-result.tau;
        double t=zriddr([&](double t){return t==right.t ? dtau-needed : integrate_panel({t,false})-needed;},
                        left.t,right.t,1e-11);
        result.tau=o.tau_target;result.length=s.x_old+t*h;result.state=path(t);result.outcome=3;
        result.quadrature_evaluations+=panels.evaluations;return finish();
      }
      result.tau+=dtau;
    }
    result.caustics+=panels.caustics;result.double_branches+=panels.double_branches;
    result.quadrature_evaluations+=panels.evaluations;
    result.turns+=(s.y_old[2]>pi/2 && s.y_new[2]<=pi/2);
    result.length=s.x_old+h;result.state=s.y_new;
    if(event>=0){result.outcome=event;break;}
    if(result.length>=o.path_limit){result.outcome=2;break;}
    s.update_old();
  }
  return finish();
}
}
