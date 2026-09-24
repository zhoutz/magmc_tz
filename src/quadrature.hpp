#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Open Gauss-Legendre rules: endpoints (including integrable singularities)
// are never sampled. The 8/16 difference is an empirical error estimator.
namespace quadrature {
template<class F> double gauss(F const &f, double a, double b, bool high, long &evals) {
  static constexpr double x8[] = {.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
  static constexpr double w8[] = {.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
  static constexpr double x16[] = {.095012509837637440185,.28160355077925891323,.45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,.94457502307323257608,.98940093499164993260};
  static constexpr double w16[] = {.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,.12462897125553387205,.095158511682492784810,.062253523938647892863,.027152459411754094852};
  const double *x = high ? x16 : x8, *w = high ? w16 : w8;
  int n = high ? 8 : 4;
  double m = (a+b)/2, h = (b-a)/2, sum = 0;
  for (int i=0; i<n; ++i) sum += w[i]*(f(m-h*x[i])+f(m+h*x[i]));
  evals += 2*n;
  return h*sum;
}
template<class F> double adaptive(F const &f, double a, double b, double atol,
                                  double rtol, long &evals, int depth=0) {
  if (a==b) return 0;
  double hi = gauss(f,a,b,true,evals), lo = gauss(f,a,b,false,evals);
  if (!std::isfinite(hi) || !std::isfinite(lo)) throw std::runtime_error("Nonfinite quadrature");
  if (std::abs(hi-lo) <= atol+rtol*std::abs(hi)) return hi;
  if (depth>=40) throw std::runtime_error("Quadrature did not converge");
  double m=(a+b)/2;
  return adaptive(f,a,m,atol/2,rtol,evals,depth+1)+adaptive(f,m,b,atol/2,rtol,evals,depth+1);
}
}
