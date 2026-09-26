// Ported from 9204013bc5919c2edad4e53c3c58bfe8cde0a4db.
#include "benchmark_io.hpp"
#include <format>
#include <map>

// Historical bench_distribution.hpp
#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../photon.hpp"
#include "../roots.hpp"
#include "../qags.hpp"
#include <array>
#include <functional>
#include <string>

inline constexpr double R_star=10, rs=1.4*schwarzschild_radius_of_sun_in_km;
inline constexpr double B_pole=1e14;

namespace transport {
// General velocity PDF interface. The solver never inspects a temperature or
// assumes a functional form. Supply ALL discontinuities, support edges and
// landmarks resolving narrow peaks in knots. A finite sampler cannot discover
// an arbitrarily narrow, unadvertised peak in a black-box function.
// bands are {low, median/peak, high} for FT07's controller, one per peak.
struct Distribution {
  std::string name;
  std::function<double(double)> pdf;
  double low, high, mean;
  std::vector<double> knots;
  std::vector<std::array<double,3>> bands;
  std::vector<double> reference_knots; // optional coarser partition for raw QAGS
  double f(double beta) const {
    return beta>low && beta<high ? pdf(beta) : 0.;
  }
  double b_bar() const { return mean; }
};

// Only this benchmark adapter knows about Boltzmann. Quantiles are cached by
// the driver, outside the timed transport loops. No distribution is truncated.
inline Distribution boltzmann_distribution(double b0) {
  Boltzmann fb(b0, 4), positive(std::abs(b0), 4);
  auto quantile=[&](double probability) {
    return zriddr([&](double b) {
      return qags([&](double v){return positive.f(v);},0.,b,1e-13,1e-12)-probability;
    },0.,1.,1e-14);
  };
  double a=quantile(.001), m=quantile(.5), b=quantile(.999);
  std::vector<double> knots={0};
  // The PDF supplies these landmarks; they are not part of the algorithm.
  for(double prob:{.001,.01,.1,.25,.5,.75,.9,.99,.999})
    knots.push_back(std::copysign(quantile(prob),b0));
  for(double fraction:{.1,.25,.5,.75,.9,.99,.999})
    knots.push_back(std::copysign(b+(1-b)*fraction,b0));
  std::array<double,3> band=b0<0 ? std::array<double,3>{-b,-m,-a} : std::array<double,3>{a,m,b};
  return {"boltzmann",[fb](double v){return fb.f(v);},fb.b_min,fb.b_max,fb.b_mean,knots,{band},{0}};
}

inline Distribution hat_distribution(bool bimodal=false) {
  // A 2e-4-wide peak, or two disconnected peaks with 20% and 80% weights.
  std::vector<std::array<double,3>> parts=bimodal ?
    std::vector<std::array<double,3>>{{-.701,-.699,.2},{-.1201,-.1199,.8}} :
    std::vector<std::array<double,3>>{{-.2001,-.1999,1.}};
  Distribution d{bimodal ? "two_hats" : "narrow_hat",{},-1,0,0,{}, {},{}};
  for(auto [lo,hi,w]:parts) {
    d.mean+=w*(lo+hi)/2;
    d.knots.insert(d.knots.end(),{lo,(lo+hi)/2,hi});
    d.bands.push_back({lo+.001*(hi-lo),(lo+hi)/2,hi-.001*(hi-lo)});
  }
  d.pdf=[parts](double v){double f=0;for(auto [lo,hi,w]:parts)if(v>lo && v<hi)f+=w/(hi-lo);return f;};
  return d;
}

inline Distribution compact_distribution() {
  // Smooth compact, asymmetric mixture on BOTH signs of beta; each component
  // is 15/(16*w)*(1-z^2)^2 for |z|<1. Integral=1, mean=center.
  const std::vector<std::array<double,3>> parts={{-.35,.0002,.8},{.5,.03,.2}};
  Distribution d{"compact_mixture",{},-1,1,-.18,{}, {},{}};
  for(auto [c,w,weight]:parts) {
    for(double z:{-1.,-.5,0.,.5,1.})d.knots.push_back(c+w*z);
    d.bands.push_back({c-w,c,c+w});
  }
  d.pdf=[parts](double v){double f=0;for(auto [c,w,weight]:parts){double z=(v-c)/w;if(std::abs(z)<1)f+=weight*15/(16*w)*std::pow(1-z*z,2);}return f;};
  return d;
}
}

// Historical bench_physics.hpp
// Adapted from origin/codex-run2; benchmark-local implementation.
#include <algorithm>

namespace transport {
using State = std::array<double,3>;

inline State geodesic(State const &y) {
  double r=y[0], alpha=y[2], L=std::sqrt(1-rs/r);
  return {L*std::cos(alpha),std::sin(alpha)/r,
          -std::sin(alpha)/(r*L)*(1-1.5*rs/r)};
}

struct ResonancePoint {
  double logx, mu, prefactor;
  double discriminant() const { return std::expm1(2*logx)+mu*mu; }
  double surface(double beta) const {
    return logx-std::log1p(-beta*mu)+.5*std::log1p(-beta*beta);
  }
};

struct Medium {
  BField const &field;
  Distribution const &distribution;
  Photon const &photon;
  mutable long evaluations=0;
  mutable long opacity_evaluations=0;

  ResonancePoint at(State const &y) const {
    ++evaluations;
    auto [r,psi,alpha]=y;
    double3 rhat=std::cos(psi)*photon.e1+std::sin(psi)*photon.e2;
    // Unit-vector roundoff must not place a pole outside the magnetic table.
    double muz=std::clamp(rhat.z,-1.,1.);
    auto B=field.calc_B(r,muz);double norm=B.length();auto b=B/norm;
    double rho=std::hypot(rhat.x,rhat.y);
    double3 theta=rho>1e-14 ? double3{rhat.x*muz/rho,rhat.y*muz/rho,-rho} : double3{std::copysign(1.,muz),0,0};
    double3 phi=rho>1e-14 ? double3{-rhat.y/rho,rhat.x/rho,0} : double3{0,1,0};
    double mu=b.x*std::cos(alpha)+std::sin(alpha)*(b.y*dot(photon.n,phi)-b.z*dot(photon.n,theta));
    return {std::log(norm*B_to_omega/photon.omega_inf)+.5*std::log1p(-rs/r),
            std::clamp(mu,-1.,1.),
            (field.p+1)*pi*field.Bphi_over_Btheta(muz)/(r*std::abs(distribution.b_bar()))};
  }

  std::array<double,2> branches(ResonancePoint const &p,double D) const {
    if(!(D>0))return {NAN,NAN};
    double x=std::exp(p.logx), den=1+D, s=x*std::sqrt(D);
    // Stable quadratic roots: the smaller numerator is recovered by Vieta.
    double large=(p.mu+std::copysign(s,p.mu))/den;
    // Use the same D as the caustic regularization, including in Vieta's
    // product. Otherwise one root retains the cancelled, unregularized D.
    double small=(p.mu*p.mu-D)/(den*large);
    if(!std::isfinite(small))small=(p.mu-std::copysign(s,p.mu))/den;
    return {std::min(large,small),std::max(large,small)};
  }

  double opacity(ResonancePoint const &p,double D,double low=-1,double high=1) const {
    ++opacity_evaluations;
    if(!(D>0) || p.prefactor==0)return 0;
    double x=std::exp(p.logx),s=std::sqrt(D), sum=0;
    // |mu-beta| = sqrt(D)*(1-beta*mu)/x. This form is stable
    // when both roots approach mu; both contributing branches are retained.
    double overlap=photon.pol==Polarization::E ? .5 : .5*D/(x*x);
    for(double b:branches(p,D)) {
      if(b<low || b>high)continue;
      double f=distribution.f(b);
      if(f==0)continue;
      sum+=f*overlap*(1-b*p.mu)*(1-b*b)*x/s;
    }
    double value=p.prefactor*sum;
    if(!std::isfinite(value) || value<0)throw std::runtime_error("Invalid resonant opacity");
    return value;
  }

  // Directional derivatives along the full geodesic, including variation of
  // B, mu, photon direction and gravitational frequency. No radial formula.
  std::array<double,2> rates(State const &y) const {
    auto dy=geodesic(y);double h=1e-5*y[0];State a=y,b=y;
    for(int i=0;i<3;++i){a[i]-=h*dy[i];b[i]+=h*dy[i];}
    auto pa=at(a),pb=at(b);
    return {(pb.logx-pa.logx)/(2*h),(pb.mu-pa.mu)/(2*h)};
  }
};

inline Photon make_photon(double r,double muz,double alpha,double azimuth,double energy,Polarization pol) {
  double sth=std::sqrt(std::max(0.,1-muz*muz));
  double3 er{sth,0,muz},etheta{muz,0,-sth},ephi{0,1,0};
  double3 tangent=std::cos(azimuth)*etheta+std::sin(azimuth)*ephi;
  return {cross(er,tangent),er,tangent,r,0,alpha,energy,pol};
}
}

// Historical bench_quadrature.hpp
// Adapted from origin/codex-run2; benchmark-local implementation.
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Open Gauss-Legendre rules: endpoints (including integrable singularities)
// are never sampled. The 8/16 difference is an empirical error estimator.
namespace quadrature {
template<class F> double gauss(F const &f, double a, double b, bool high, long &evals) {
  static constexpr double x8[] = {.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
  static constexpr double w8[] = {.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
  static constexpr double x16[] = {.095012509837637440185,.28160355077925891323,.45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,.94457502307323257608,.98940093499164993260};
  static constexpr double w16[] = {.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,.12462897125553387205,.095158511682492784810,.062253523938647892863,.027152459411754094852};
  const double *x = high ? x16 : x8, *w = high ? w16 : w8;
  int n = high ? 8 : 4;
  double m = (a+b)/2, h = (b-a)/2, sum = 0;
  for (int i=0; i<n; ++i) sum += w[i]*(f(m-h*x[i])+f(m+h*x[i]));
  evals += 2*n;
  return h*sum;
}
template<class F> double adaptive(F const &f, double a, double b, double atol,
                                  double rtol, long &evals, int depth=0) {
  if (a==b) return 0;
  double hi = gauss(f,a,b,true,evals), lo = gauss(f,a,b,false,evals);
  if (!std::isfinite(hi) || !std::isfinite(lo)) throw std::runtime_error("Nonfinite quadrature");
  if (std::abs(hi-lo) <= atol+rtol*std::abs(hi)) return hi;
  if (depth>=40) throw std::runtime_error("Quadrature did not converge");
  double m=(a+b)/2;
  return adaptive(f,a,m,atol/2,rtol,evals,depth+1)+adaptive(f,m,b,atol/2,rtol,evals,depth+1);
}
}

// Historical bench_panels.hpp
// Adapted from origin/codex-run2; benchmark-local implementation.
#include "../roots.hpp"
#include <vector>

namespace transport {
struct Boundary { double t; bool singular=false; };

// A path segment is parameterized by t in [0,1], NOT by radius. Events
// remain valid through radial turning points and changing magnetic angles.
template<class Path, class Material=Medium> struct ResonancePanels {
  Material const &medium;
  Path const &path;
  double length, low, high;
  std::vector<Boundary> boundaries;
  long evaluations=0;
  int caustics=0, double_branches=0;

  ResonancePoint point(double t) const { return medium.at(path(t)); }
  double surface(double t,double beta) const {
    auto p=point(t);
    return beta==2 ? p.logx-.5*std::log(std::max(1e-30,1-p.mu*p.mu)) : p.surface(beta);
  }

  ResonancePanels(Material const &m,Path const &p,double h,std::vector<double> const &velocities,
                  int probes,double lo=-1,double hi=1)
      :medium(m),path(p),length(h),low(lo),high(hi) {
    boundaries={{0,false},{1,false}};
    std::vector<ResonancePoint> samples;
    for(int i=0;i<=probes;++i)samples.push_back(point(double(i)/probes));
    std::vector<double> surfaces=velocities;surfaces.push_back(2); // D=0
    for(double beta:surfaces) {
      std::vector<std::pair<double,double>> nodes;
      for(int i=0;i<=probes;++i) {
        auto p=samples[i];
        nodes.push_back({double(i)/probes,beta==2 ? p.logx-.5*std::log(std::max(1e-30,1-p.mu*p.mu)) : p.surface(beta)});
      }
      // Endpoints can have equal signs around a thin resonance pocket.
      // Locate interior extrema before root bracketing, in addition to
      // sampling. Probe refinement is independently checked in validation.
      for(int i=1;i<probes;++i) {
        double s1=nodes[i].second-nodes[i-1].second,s2=nodes[i+1].second-nodes[i].second;
        if(s1*s2>=0)continue;
        double a=nodes[i-1].first,b=nodes[i+1].first,sgn=s1<0 ? 1. : -1.;
        constexpr double g=.6180339887498948482;
        double c=b-g*(b-a),d=a+g*(b-a),fc=sgn*surface(c,beta),fd=sgn*surface(d,beta);
        for(int k=0;k<48;++k) {
          if(fc<fd){b=d;d=c;fd=fc;c=b-g*(b-a);fc=sgn*surface(c,beta);}
          else {a=c;c=d;fc=fd;d=a+g*(b-a);fd=sgn*surface(d,beta);}
        }
        double t=(a+b)/2;
        nodes.push_back({t,surface(t,beta)});
      }
      std::sort(nodes.begin(),nodes.end());
      for(size_t i=1;i<nodes.size();++i) {
        auto [a,fa]=nodes[i-1];auto [b,fb]=nodes[i];
        if(a==b)continue;
        if(fa==0)boundaries.push_back({a,beta==2});
        if(fb==0)boundaries.push_back({b,beta==2});
        if(std::signbit(fa)==std::signbit(fb))continue;
        double t=zriddr([&](double t){return surface(t,beta);},a,b,3e-15);
        boundaries.push_back({t,beta==2});
      }
    }
    std::sort(boundaries.begin(),boundaries.end(),[](auto a,auto b){return a.t<b.t;});
    std::vector<Boundary> clean;
    for(auto b:boundaries) {
      if(!clean.empty() && b.t-clean.back().t<8e-15)clean.back().singular|=b.singular;
      else clean.push_back(b);
    }
    boundaries=std::move(clean);
    for(auto b:boundaries)caustics+=b.singular;
  }

  double integrate(Boundary left,Boundary right,double atol,double rtol) {
    double a=left.t,b=right.t,w=b-a;
    if(!(w>0))return 0;
    auto mid=point((a+b)/2);double D=mid.discriminant();
    if(D<=0)return 0; // all D=0 crossings have been split
    int contributing=0;
    for(double beta:medium.branches(mid,D))contributing+=beta>std::max(low,medium.distribution.low) && beta<std::min(high,medium.distribution.high);
    if(contributing==0)return 0;
    double_branches+=contributing==2;

    // Local Taylor coefficients anchor D at a located simple caustic.
    // They avoid subtracting nearly equal O(1) numbers in tiny neighborhoods
    // after the square-root change of variable. No angular specialization.
    auto coefficients=[&](double t,double direction) {
      auto y=path(t),p=medium.at(y);auto rates=medium.rates(y);
      double first=2*std::exp(2*p.logx)*rates[0]+2*p.mu*rates[1];
      double dl=1e-4*y[0];auto dy=geodesic(y);State yp=y,ym=y;
      for(int k=0;k<3;++k){yp[k]+=dl*dy[k];ym[k]-=dl*dy[k];}
      auto pp=medium.at(yp),pm=medium.at(ym);
      auto rp=medium.rates(yp),rm=medium.rates(ym);
      double second=(2*std::exp(2*pp.logx)*rp[0]+2*pp.mu*rp[1]
                    -2*std::exp(2*pm.logx)*rm[0]-2*pm.mu*rm[1])/(2*dl);
      return std::array<double,2>{direction*first*length,second*length*length};
    };
    std::array<double,2> lc{},rc{};
    if(left.singular)lc=coefficients(a,1);
    if(right.singular)rc=coefficients(b,-1);
    auto integrand=[&](double z) {
      // sin^2 maps both possible endpoint singularities at once.
      double sa=std::sin(.5*pi*z),sb=std::cos(.5*pi*z);
      double da=w*sa*sa,db=w*sb*sb;
      double t=z<.5 ? a+da : b-db;
      auto p=point(t);double D=p.discriminant();
      // A panel may be extremely narrow in D. Scaling the Taylor region by
      // panel width alone leaves most of such a panel dominated by roundoff.
      if(left.singular && da*length<1e-5*path(a)[0] && std::abs(lc[0]*da)<1e-5)
        D=lc[0]*da+.5*lc[1]*da*da;
      if(right.singular && db*length<1e-5*path(b)[0] && std::abs(rc[0]*db)<1e-5)
        D=rc[0]*db+.5*rc[1]*db*db;
      return medium.opacity(p,D,low,high)*length*w*pi*sa*sb;
    };
    return quadrature::adaptive(integrand,0,1,atol,rtol,evaluations);
  }
};
}

// Historical bench_transport.hpp
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

namespace velocity_panels_bench {
inline BField field("table/bfield_t10.txt",B_pole,R_star);
inline std::map<double,transport::Distribution> cache;
void prepare(std::vector<optical_bench::Input> const &cases,std::string const &mode) {
  if(mode=="cached")
    for(auto const &c:cases)
      if(!cache.contains(c.b0)) cache.emplace(c.b0,transport::boltzmann_distribution(c.b0));
}
double total_optical_depth(optical_bench::Input const &c,std::string const &mode,
                           optical_bench::Counts &work) {
  if(mode!="default" && mode!="matched" && mode!="cached" && mode!="audit")
    throw std::invalid_argument("Unknown velocity_panels configuration");
  auto start=optical_bench::Clock::now();
  // Default includes all distribution preparation, like the other candidates.
  auto distribution=mode=="cached" ? cache.at(c.b0) : transport::boltzmann_distribution(c.b0);
  work.preparation_seconds=std::chrono::duration<double>(optical_bench::Clock::now()-start).count();
  auto photon=transport::make_photon(c.radius,c.muz,c.alpha,c.az,c.energy,c.polarization());
  transport::Options cfg;
  if(mode=="matched") cfg.tolerance=1e-8;
  if(mode=="audit") {cfg.tolerance=1e-9;cfg.resolution=.05;cfg.probes=16;}
  auto r=transport::integrate(field,distribution,photon,transport::Method::panels,cfg);
  work.geometry=r.evaluations;work.density=r.opacity_evaluations;work.steps=r.steps;
  // The historical implementation counts nodes, not quadrature invocations.
  work.quadrature=r.quadrature_evaluations;
  return r.tau;
}
}
#ifndef OPTICAL_BENCH_NO_MAIN
int main(int argc,char **argv) {
  return optical_bench::run(argc,argv,"velocity_panels",velocity_panels_bench::total_optical_depth,
                            velocity_panels_bench::prepare);
}
#endif
