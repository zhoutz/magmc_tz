#pragma once
#include "../src/transport.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
struct Case { double b0,muz; int pol; double energy,reference; };
struct Result { double tau=0; long evaluations=0,steps=0,rejected=0; std::string status="ok"; };
struct Config { double tol=1e-6,parameter=.1; int repeats=5; double initial_step=.01; };
using Edges=std::array<double,3>;
inline Edges distribution_edges(Boltzmann const &fb){return transport::distribution_edges(fb);}
inline Result general_solve(BField const &field,Case const &c,Config const &cfg,Edges const &edges,transport::Method method) {
 Boltzmann fb(c.b0);
 auto photon=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
 transport::Options options;options.method=method;options.tolerance=cfg.tol;options.resolution=cfg.parameter;options.initial_step=cfg.initial_step;
 auto r=transport::integrate(field,fb,photon,options,&edges);
 return {r.tau,r.evaluations,r.steps,r.quadrature_evaluations,"ok"};
}
inline Result baseline(BField const &field,Case const &c,Config const &cfg) {
 Boltzmann fb(c.b0);auto p=transport::make_photon(R_star,c.muz,0,0,c.energy,c.pol ? Polarization::E : Polarization::O);
 PhotonEvolution ev{field,fb,p.n,p.e1,p.e2,{R_star,0,0,0},c.energy,p.pol};
 long calls=0;auto rhs=[&](double l,YVector const &y,YVector &dy){++calls;ev(l,y,dy);};
 StepperDopr5<4,decltype(rhs)> s(rhs,cfg.tol,cfg.tol);s.add_event([](double,YVector const &y){return y[0]-10000;});s.init(0,cfg.initial_step,{R_star,0,0,0});
 long steps=0;while(true){if(++steps>2000000)throw std::runtime_error("Step budget exceeded");s.do_step();if(s.detect_event()==0)break;s.update_old();}
 return {s.y_new[3],calls,steps,0,"ok"};
}
template<class Solver> int benchmark_main(int argc,char **argv,std::string name,Solver solve) {
 try {
  if (argc<2) throw std::runtime_error("Usage: executable output.csv [tolerance] [parameter] [repeats] [initial_step]");
  Config cfg;
  if(argc>2) cfg.tol=std::stod(argv[2]);
  if(argc>3) cfg.parameter=std::stod(argv[3]);
  if(argc>4) cfg.repeats=std::stoi(argv[4]);
  if(argc>5) cfg.initial_step=std::stod(argv[5]);
  if (!(cfg.tol>0 && cfg.parameter>0 && cfg.repeats>0 && cfg.initial_step>0)) throw std::runtime_error("Invalid config");
  BField field("table/bfield_t10.txt",B_pole,R_star);
  std::ifstream input("table/bench_od.txt");
  if(!input) throw std::runtime_error("Missing reference file");
  std::vector<Case> cases; std::string line;
  while(std::getline(input,line)) { if(line.empty()||line[0]=='#')continue; Case c; std::istringstream in(line); if(!(in>>c.b0>>c.muz>>c.pol>>c.energy>>c.reference))throw std::runtime_error("Malformed reference"); cases.push_back(c); }
  if(cases.size()!=900) throw std::runtime_error("Expected 900 reference cases");
  auto setup_start=std::chrono::steady_clock::now();
  std::map<double,Edges> edges;
  for(auto const &c:cases) if(!edges.count(c.b0)) edges[c.b0]=distribution_edges(Boltzmann(c.b0));
  double setup_ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-setup_start).count();
  auto safe_solve=[&](Case const &c) {
    try { return solve(field,c,cfg,edges.at(c.b0)); }
    catch(std::exception const &e) { return Result{std::numeric_limits<double>::quiet_NaN(),0,0,0,e.what()}; }
  };
  // Warm up before timing; every timed batch performs all 900 solves, no file I/O.
  std::vector<Result> results;
  for(auto const &c:cases) results.push_back(safe_solve(c));
  std::vector<std::vector<double>> elapsed(cases.size());
  std::vector<double> batches;
  for(int rep=0;rep<cfg.repeats;++rep) {
    auto batch_start=std::chrono::steady_clock::now();
    for(size_t i=0;i<cases.size();++i) {
      auto start=std::chrono::steady_clock::now();
      Result r=safe_solve(cases[i]);
      elapsed[i].push_back(std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-start).count());
      if(r.status!=results[i].status || (r.status=="ok" && r.tau!=results[i].tau)) throw std::runtime_error("Nondeterministic result");
    }
    batches.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-batch_start).count());
  }
  std::sort(batches.begin(),batches.end());
  std::ofstream out(argv[1]); if(!out)throw std::runtime_error("Cannot open output");
  out<<std::setprecision(17);
  out<<"# method="<<name<<" tolerance="<<cfg.tol<<" parameter="<<cfg.parameter<<" initial_step="<<cfg.initial_step<<" repeats="<<cfg.repeats<<" setup_ms="<<setup_ms<<" batch_median_ms="<<batches[batches.size()/2]<<" batch_min_ms="<<batches.front()<<" batch_max_ms="<<batches.back()<<'\n';
  out<<"b0,muz,pol,omega_inf,tau,reference,abs_error,rel_error,microseconds,evaluations,steps,quadrature_evaluations,status\n";
  double maxrel=0; int failed=0, exceptions=0;
  for(size_t i=0;i<cases.size();++i) {
    auto c=cases[i]; auto r=results[i]; auto &times=elapsed[i]; std::sort(times.begin(),times.end());
    double ae=std::abs(r.tau-c.reference), re=ae/std::abs(c.reference);
    maxrel=std::max(maxrel,re); failed+=re>1e-3;
    exceptions+=r.status!="ok";
    out<<c.b0<<','<<c.muz<<','<<c.pol<<','<<c.energy<<','<<r.tau<<','<<c.reference<<','<<ae<<','<<re<<','<<times[times.size()/2]<<','<<r.evaluations<<','<<r.steps<<','<<r.rejected<<','<<r.status<<'\n';
  }
  std::cout<<name<<" tol="<<cfg.tol<<" param="<<cfg.parameter<<" maxrel="<<maxrel<<" errors>0.1%="<<failed<<" failures="<<exceptions<<" batch_ms="<<batches[batches.size()/2]<<'\n';
  return 0;
 } catch(std::exception const &e) { std::cerr<<name<<": "<<e.what()<<'\n'; return 1; }
}
