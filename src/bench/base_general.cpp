// General-direction version of the unguarded coupled RK method in base.cpp.
#include "bench_driver.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"base_general",transport::Method::baseline);}
