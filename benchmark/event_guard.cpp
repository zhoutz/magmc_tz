#include "benchmark_common.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"event_guard",[](auto const &f,auto const &c,auto const &cfg,auto const &e){return general_solve(f,c,cfg,e,transport::Method::event_guard);});}
