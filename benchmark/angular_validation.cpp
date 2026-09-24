#include "../output/benchmark_common.hpp"
#include <random>

int main(int argc,char **argv) {
 int n=argc>1 ? std::stoi(argv[1]) : 80;
 std::ofstream out("output/angular_validation.csv");out<<std::setprecision(17);
 out<<"id,b0,muz,alpha,azimuth,r,energy,pol,method,tau,reference,relative_error,ms,steps,caustics,double_branches,turns,outcome,status\n";
 BField field("table/bfield_t10.txt",B_pole,R_star);
 std::mt19937_64 rng(20260924);std::uniform_real_distribution<double> U(0,1);
 int failures=0;
 for(int i=0;i<n;++i) {
  double b0=std::array<double,4>{-.1,-.5,-.9,.3}[i%4];
  double muz=2*U(rng)-1,alpha=pi*U(rng),az=2*pi*U(rng),r=30,energy=std::pow(10.,-1+3*U(rng));
  // Explicit pole, equator, radial, tangent, inward and oblique cases.
  if(i<16){muz=std::array<double,4>{-1.,-.5,0.,1.}[i%4];alpha=std::array<double,4>{0.,pi/2,2.7,3.05}[i/4];az=.7;energy=1;}
  if(i>=72){r=R_star;alpha=std::array<double,4>{.2,.6,pi/2,3.05}[i%4];}
  Boltzmann fb(b0);auto edges=transport::distribution_edges(fb);
  for(auto pol:{Polarization::E,Polarization::O}) {
   auto photon=transport::make_photon(r,muz,alpha,az,energy,pol);
   transport::Options ref;ref.method=transport::Method::event_mesh;ref.resolution=.005;ref.tolerance=1e-8;ref.probes=16;
   transport::Result truth;
   try {truth=transport::integrate(field,fb,photon,ref);}
   catch(std::exception const &e){std::cerr<<"reference id="<<i<<" pol="<<int(pol)<<": "<<e.what()<<'\n';++failures;continue;}
   for(auto method:{transport::Method::event_guard,transport::Method::phase_cap,transport::Method::ft07}) {
    std::string name=method==transport::Method::ft07 ? "ft07" : method==transport::Method::phase_cap ? "phase_cap" : method==transport::Method::event_mesh ? "event_mesh" : "event_guard";
    transport::Options cfg;cfg.method=method;cfg.tolerance=1e-6;cfg.resolution=method==transport::Method::ft07 ? .01 : .1;
    auto start=std::chrono::steady_clock::now();
    try {
     auto result=transport::integrate(field,fb,photon,cfg,&edges);
     double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
     double reference=truth.tau;
     if(method==transport::Method::ft07){auto truncated=ref;truncated.truncate_distribution=true;reference=transport::integrate(field,fb,photon,truncated,&edges).tau;}
     double error=std::abs(result.tau-reference)/std::max(1e-12,std::abs(reference));
     std::string status=std::abs(result.tau-reference)<=1e-7+1e-4*std::abs(reference) && result.outcome==truth.outcome ? "ok" : "mismatch";
     failures+=status!="ok";
     out<<i<<','<<b0<<','<<muz<<','<<alpha<<','<<az<<','<<r<<','<<energy<<','<<int(pol==Polarization::E)<<','<<name<<','<<result.tau<<','<<reference<<','<<error<<','<<ms<<','<<result.steps<<','<<result.caustics<<','<<result.double_branches<<','<<result.turns<<','<<int(result.outcome)<<','<<status<<'\n';
     if(status!="ok")std::cerr<<"mismatch id="<<i<<" mode="<<name<<" pol="<<int(pol)<<" tau="<<result.tau<<" ref="<<reference<<" error="<<error<<'\n';
    }catch(std::exception const &e){++failures;out<<i<<','<<b0<<','<<muz<<','<<alpha<<','<<az<<','<<r<<','<<energy<<','<<int(pol==Polarization::E)<<','<<name<<",nan,"<<truth.tau<<",nan,0,0,0,0,0,0,"<<e.what()<<'\n';std::cerr<<"id="<<i<<" mode="<<name<<" pol="<<int(pol)<<": "<<e.what()<<'\n';}
    out.flush();
   }
  }
  std::cout<<"angular case "<<i+1<<'/'<<n<<" failures="<<failures<<std::endl;
 }
 return failures ? 1 : 0;
}
