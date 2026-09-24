#include "../src/transport_fast.hpp"
#include "../src/resonance.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
int main(int argc,char **argv) {
 try {
 if(argc!=5)throw std::runtime_error("Usage: executable radial|angular id rtol max_intervals");
 std::string suite=argv[1];int wanted_id=std::stoi(argv[2]);
 fast_transport::TransportOptions o;o.guard_extremum_bounds=true;o.quadrature_rtol=std::stod(argv[3]);o.quadrature_atol=.01*o.quadrature_rtol;o.max_quadrature_intervals=std::stoull(argv[4]);
 std::ifstream input("output/quad_refactor_cases.txt");std::string line;
 while(std::getline(input,line)) {
  if(line.empty()||line[0]=='#')continue;
  std::istringstream row(line);std::string group;int id,pol,outcome;double b0,muz,alpha,azimuth,r,psi,energy,truth,initial,target;
  row>>group>>id>>b0>>muz>>alpha>>azimuth>>r>>psi>>energy>>pol>>truth>>outcome>>initial>>target;
  if(group!=suite||id!=wanted_id)continue;
  BField field("table/bfield_t10.txt",B_pole,R_star);Boltzmann fb(b0);auto p=transport::make_photon(r,muz,alpha,azimuth,energy,Polarization(pol));p.psi=psi;
  double tau=NAN;std::string status="ok";auto start=std::chrono::steady_clock::now();
  try{tau=fast_transport::propagate(field,fb,p,o).tau;}catch(std::exception const &e){status=e.what();}
  double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
  std::cout<<std::setprecision(17)<<suite<<','<<id<<','<<o.quadrature_rtol<<','<<o.max_quadrature_intervals<<','<<tau<<','<<truth<<','<<ms<<','<<std::quoted(status)<<'\n';return 0;
 }
 throw std::runtime_error("Unknown case");
 }catch(std::exception const &e){std::cerr<<e.what()<<'\n';return 1;}
}
