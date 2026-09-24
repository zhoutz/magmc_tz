// Time the UNMODIFIED total_optical_depth function, excluding startup and I/O.
// Compile with -DTIME_REF for ref.cpp, otherwise use base.cpp.
#define main original_main
#ifdef TIME_REF
#include "ref.cpp"
#else
#include "base.cpp"
#endif
#undef main
#include <chrono>
#include <fstream>
#include <iomanip>
#include <vector>
int main() {
#ifdef TIME_REF
  const char *name="ref";
#else
  const char *name="base";
#endif
  struct Case {double b0,muz,energy;int pol;double truth;};
  std::ifstream in("table/bench_od.txt");std::vector<Case> cases;Case c;
  while(in>>c.b0>>c.muz>>c.energy>>c.pol>>c.truth)cases.push_back(c);
  std::vector<double> values,times;
  for(int rep=0;rep<4;++rep) {
    auto t=std::chrono::steady_clock::now();
    values.clear();
    for(auto c:cases)values.push_back(total_optical_depth(c.b0,c.muz,c.energy,c.pol ? Polarization::E : Polarization::O));
    if(rep)times.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count());
  }
  std::sort(times.begin(),times.end());
  std::ofstream data(std::string("output/")+name+"_original.txt"),stats(std::string("output/")+name+"_original_stats.txt");
  data<<std::setprecision(17);stats<<std::setprecision(17);
  stats<<"# method="<<name<<" suite=radial repeats=3 median_ms="<<times[1]<<" min_ms="<<times[0]<<" max_ms="<<times[2]<<'\n';
  for(size_t i=0;i<cases.size();++i){auto c=cases[i];data<<c.b0<<' '<<c.muz<<' '<<c.energy<<' '<<c.pol<<' '<<values[i]<<'\n';}
  std::cout<<name<<" original median_ms="<<times[1]<<'\n';
}
