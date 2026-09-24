// Generalization of ref.cpp: small orbit steps + spatial QAGS, including all
// support boundaries and continuing after a D=0 crossing until escape/star.
#include "bench_driver.hpp"
int main(int argc,char **argv){return benchmark_main(argc,argv,"ref_general",transport::Method::reference);}
