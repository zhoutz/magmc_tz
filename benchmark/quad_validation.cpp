#include "../src/quad.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>

int main() {
 try {
  std::ofstream out("output/quad_validation.csv");out<<std::setprecision(17)<<"case,integral,exact,absolute_error,estimated_error,evaluations,panels\n";
  int tests=0;
  auto check=[&](char const *name,auto const &f,double a,double b,double exact,double tolerance,quad::Options options=quad::Options{}) {
   std::size_t calls=0;quad::Statistics stats;
   auto q=quad::integrate([&](double x){++calls;return f(x);},a,b,options,stats);
   out<<name<<','<<q.integral<<','<<exact<<','<<std::abs(q.integral-exact)<<','<<q.error<<','<<stats.evaluations<<','<<stats.intervals<<'\n';
   if(std::abs(q.integral-exact)>tolerance || calls!=stats.evaluations || q.error<0)throw std::runtime_error(std::string("Quadrature mismatch: ")+name);
   ++tests;
  };
  check("polynomial",[](double x){return x*x*x*x;},0,2,6.4,1e-13);
  check("signed",[](double x){return x*x-1;},-1,1,-4./3,1e-13);
  check("reversed",[](double x){return std::sin(x);},std::numbers::pi,0,-2,1e-13);
  check("empty",[](double){throw std::runtime_error("empty interval evaluated");return 0.;},1,1,0,0);
  check("left_sqrt_singularity",[](double x){return 1/std::sqrt(x);},0,1,2,5e-7);
  check("right_sqrt_singularity",[](double x){return 1/std::sqrt(1-x);},0,1,2,5e-7);
  check("both_sqrt_singularities",[](double x){return 1/std::sqrt(x*(1-x));},0,1,std::numbers::pi,8e-7);
  check("log_singularity",[](double x){return std::log(x);},0,1,-1,3e-8);
  auto expect_throw=[&](auto const &fn){bool threw=false;try{fn();}catch(std::exception const &){threw=true;}if(!threw)throw std::runtime_error("Expected quadrature failure was not reported");++tests;};
  expect_throw([]{quad::integrate([](double){return INFINITY;},0,1);});
  expect_throw([]{quad::integrate([](double x){return std::abs(x-.3);},0,1,{1e-15,1e-15,1});});
  expect_throw([]{quad::integrate([](double x){return x;},0,1,{0,0,1});});
  expect_throw([]{quad::integrate([](double x){return x;},0,INFINITY);});
  std::cout<<"quad validation: "<<tests<<" checks passed\n";
 }catch(std::exception const &e){std::cerr<<e.what()<<'\n';return 1;}
}
