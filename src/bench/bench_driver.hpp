#pragma once
#include "bench_transport.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>

struct BenchCase {
  double b0,muz,energy,r=R_star,alpha=0,azimuth=0;
  int pol=1,distribution=0;
};
inline std::vector<BenchCase> benchmark_cases(std::string const &suite) {
  std::vector<BenchCase> cases;
  if(suite=="radial") {
    std::ifstream in("table/bench_od.txt");
    if(!in)throw std::runtime_error("Missing radial table");
    BenchCase c;double tau;
    while(in>>c.b0>>c.muz>>c.energy>>c.pol>>tau)cases.push_back(c);
    if(cases.size()!=900)throw std::runtime_error("Expected 900 radial cases");
  } else if(suite=="angular") {
    // Outgoing, nearly tangent, inward then turning, and star-hitting rays.
    // Azimuth varies, so most trajectories are not magnetic meridians.
    for(double b0:{-.1,-.3,-.7,-.9})
      for(double muz:{-.8,-.3,0.,.4,.85})
        for(int geometry=0;geometry<4;++geometry)
          for(double energy:{.1,10.,100.})
            for(int pol:{1,0}) {
              double al[]={.4,1.54,2.0,2.9};
              double az[]={.31,1.7,3.4,5.1};
              cases.push_back({b0,muz,energy,geometry<2 ? 10. : 40.,al[geometry],az[geometry],pol,0});
            }
    // Fixed-seed additional scattered rays: b0 colder than the radial table,
    // arbitrary starting radii and full angular range.
    std::mt19937_64 gen(20260925);
    auto u=[&](){return std::generate_canonical<double,53>(gen);};
    for(int i=0;i<96;++i) {
      double muz=1.96*u()-.98,energy=std::pow(10.,-2+4*u());
      double r=10+90*u(),alpha=std::acos(2*u()-1),az=2*pi*u();
      cases.push_back({i%2 ? -.03 : -.5,muz,energy,r,alpha,az,i%2,0});
    }
  } else if(suite=="general") {
    for(int dist:{1,2,3})
      for(double muz:{-.7,-.2,.2,.7})
        for(double alpha:{.6,1.5,2.0})
          for(double energy:{1.,100.})
            for(int pol:{1,0})
              cases.push_back({0,muz,energy,alpha>pi/2 ? 40. : 10.,alpha,1.1+muz,pol,dist});
  } else throw std::runtime_error("Unknown suite: "+suite);
  return cases;
}

inline int benchmark_main(int argc,char **argv,std::string name,transport::Method method) {
  using namespace transport;
  try {
    gsl_set_error_handler_off();
    std::string suite=argc>1 ? argv[1] : "radial";
    std::string stem=argc>2 ? argv[2] : name+(suite=="radial" ? "" : "_"+suite);
    Options options;
    options.resolution=method==Method::ft07 ? .01 : method==Method::geometric ? .001 : method==Method::reference ? .01 : .2;
    if(argc>3)options.tolerance=std::stod(argv[3]);
    if(argc>4)options.resolution=std::stod(argv[4]);
    int repeats=argc>5 ? std::stoi(argv[5]) : 3;
    if(argc>6)options.initial_step=std::stod(argv[6]);
    if(argc>7)options.probes=std::stoi(argv[7]);
    if(repeats<1)throw std::runtime_error("repeats must be positive");
    auto cases=benchmark_cases(suite);
    BField field("table/bfield_t10.txt",B_pole,R_star);
    auto start=std::chrono::steady_clock::now();
    std::map<double,Distribution> boltz;
    std::vector<Distribution> other={hat_distribution(),hat_distribution(true),compact_distribution()};
    for(auto c:cases)if(c.distribution==0 && !boltz.count(c.b0))boltz.emplace(c.b0,boltzmann_distribution(c.b0));
    double setup=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    std::vector<Result> results(cases.size());
    std::vector<std::string> errors(cases.size());
    auto solve=[&](size_t i) {
      auto c=cases[i];
      auto const &dist=c.distribution==0 ? boltz.at(c.b0) : other.at(c.distribution-1);
      auto p=make_photon(c.r,c.muz,c.alpha,c.azimuth,c.energy,c.pol ? Polarization::E : Polarization::O);
      try {return integrate(field,dist,p,method,options);}
      catch(std::exception const &e) {
        if(errors[i].empty())std::cerr<<name<<" "<<suite<<" case "<<i<<": "<<e.what()<<'\n';
        errors[i]=e.what();Result r;r.tau=NAN;return r;
      }
    };
    // Warm-up, then time whole batches with file I/O excluded.
    for(size_t i=0;i<cases.size();++i)results[i]=solve(i);
    std::vector<double> times;
    for(int rep=0;rep<repeats;++rep) {
      start=std::chrono::steady_clock::now();
      for(size_t i=0;i<cases.size();++i) {
        auto r=solve(i);
        if(std::isfinite(r.tau) && r.tau!=results[i].tau)throw std::runtime_error("Non-deterministic result");
      }
      times.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::ofstream data("output/"+stem+".txt"),stats("output/"+stem+"_stats.txt");
    if(!data || !stats)throw std::runtime_error("Cannot open output files");
    data<<std::setprecision(17);stats<<std::setprecision(17);
    stats<<"# method="<<name<<" suite="<<suite<<" tolerance="<<options.tolerance<<" resolution="<<options.resolution
         <<" repeats="<<repeats<<" probes="<<options.probes<<" initial_step="<<options.initial_step
         <<" setup_ms="<<setup<<" median_ms="<<times[times.size()/2]<<" min_ms="<<times.front()<<" max_ms="<<times.back()<<'\n';
    stats<<"# id tau length radius steps field_evals opacity_evals quad_evals caustics double_branch_panels turns outcome status\n";
    int failed=0;
    if(suite!="radial")data<<"# id distribution b0 r muz alpha azimuth energy pol tau\n";
    for(size_t i=0;i<cases.size();++i) {
      auto c=cases[i];auto r=results[i];
      if(suite=="radial")data<<c.b0<<' '<<c.muz<<' '<<c.energy<<' '<<c.pol<<' '<<r.tau<<'\n';
      else data<<i<<' '<<c.distribution<<' '<<c.b0<<' '<<c.r<<' '<<c.muz<<' '<<c.alpha<<' '<<c.azimuth<<' '<<c.energy<<' '<<c.pol<<' '<<r.tau<<'\n';
      failed+=!errors[i].empty();
      stats<<i<<' '<<r.tau<<' '<<r.length<<' '<<r.state[0]<<' '<<r.steps<<' '<<r.evaluations<<' '<<r.opacity_evaluations<<' '
           <<r.quadrature_evaluations<<' '<<r.caustics<<' '<<r.double_branches<<' '<<r.turns<<' '<<r.outcome<<' '
           <<(errors[i].empty() ? "ok" : errors[i])<<'\n';
    }
    std::cout<<stem<<" cases="<<cases.size()<<" median_ms="<<times[times.size()/2]<<" failures="<<failed<<'\n';
    return failed ? 2 : 0;
  } catch(std::exception const &e) {std::cerr<<name<<": "<<e.what()<<'\n';return 1;}
}
