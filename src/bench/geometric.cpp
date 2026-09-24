// The smallest change to base.cpp: cap each coupled orbit/tau RK step at eta*r.
#include "bench_driver.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"geometric",transport::Method::geometric);}
