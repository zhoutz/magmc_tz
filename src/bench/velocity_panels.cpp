// Locate crossings of the distribution's velocity landmarks and D=0 on the
// dense geodesic. Integrate each panel with a sin^2 endpoint regularization.
#include "bench_driver.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"velocity_panels",transport::Method::panels);}
