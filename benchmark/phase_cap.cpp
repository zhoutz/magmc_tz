#include "../output/benchmark_common.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"phase_cap",[](auto const &f,auto const &c,auto const &cfg,auto const &e){return general_solve(f,c,cfg,e,transport::Method::phase_cap);});}
