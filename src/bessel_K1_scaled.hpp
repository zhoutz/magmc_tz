#pragma once

// bessel_K1_scaled.hpp -- C++11, header-only, standard library only.
//
// Computes exp(x) * K_1(x), where K_1 is the modified cylindrical Bessel
// function of the second kind and order one. All arithmetic uses double.
// Intended for IEEE-754 binary64 with ordinary floating-point semantics
// (not -ffast-math). No allocation, exceptions, or mutable shared state.
//
// Usage:
//   #include "bessel_K1_scaled.hpp"
//   double scaled = bessel_K1_scaled(1000.0);
//   double ratio  = 1.0 / scaled; // exp(-1000) / K_1(1000)
//
// Domain / special values:
//   finite x > 0 : exp(x) * K_1(x)
//   x == +/-0    : +infinity (right-hand limit)
//   x < 0        : quiet NaN (including -infinity)
//   x == +inf    : +0 (limit as x -> +infinity)
//   x == NaN     : NaN
// For extremely small positive x the mathematical result can exceed
// DBL_MAX; in that case the result is +infinity. No explicit errno or
// floating-point exception handling is performed.
//
// The rational approximations / coefficients are adapted from the
// 53-bit implementation in Boost.Math 1.76.0:
// https://live.boost.org/doc/libs/1_76_0/boost/math/special_functions/detail/bessel_k1.hpp
// Changes: standalone double-only evaluation, removal of exp(-x) from
// the large-x branch, scaling of the small-x branch, a tiny-x branch,
// and explicit handling of domain / special values.
// This file does not require Boost or GSL headers or libraries.
// The original copyright notices and license are reproduced below.
//
// Validation (GCC, C++11, -O3, binary64): 380,768 sampled inputs compared
// with SciPy k1e; 1,117 inputs also checked against independent formulas
// evaluated with 100-digit arithmetic. Maximum observed relative error
// in the latter check was 3.6e-16 for finite results. This is empirical
// validation, not a proof of correct rounding or a uniform error bound.

#include <cmath>
#include <cstddef>
#include <limits>

namespace bessel_k1_scaled_detail {

// Coefficients are in ascending powers: c[0] + c[1]*z + ... .
// Every call below has 0 <= z <= 1, so ordinary Horner evaluation
// suffices; reversing the polynomials for large z is unnecessary.
template <std::size_t N>
inline double rational(const double (&p)[N], const double (&q)[N],
                       double z) noexcept {
  double num = p[N - 1];
  double den = q[N - 1];
  for (std::size_t i = N - 1; i > 0; --i) {
    num = num * z + p[i - 1];
    den = den * z + q[i - 1];
  }
  return num / den;
}

} // namespace bessel_k1_scaled_detail

inline double bessel_K1_scaled(double x) noexcept {
  if (std::isnan(x)) {
    return x;
  }
  if (x < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (x == 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  if (std::isinf(x)) {
    return 0.0;
  }

  // exp(x)*K_1(x) = 1/x + 1 + O(x*log(x)).
  // Here the omitted correction is far below double precision.
  // This also avoids underflow in the small-x polynomial arguments.
  if (x < std::numeric_limits<double>::epsilon()) {
    return 1.0 / x + 1.0;
  }

  using bessel_k1_scaled_detail::rational;

  if (x <= 1.0) {
    // Approximate I_1(x), the coefficient of log(x) in K_1(x).
    static const double p_i[] = {
        -3.62137953440350228e-03, 7.11842087490330300e-03,
        1.00302560256614306e-05, 1.77231085381040811e-06};
    static const double q_i[] = {
        1.00000000000000000e+00, -4.80414794429043831e-02,
        9.85972641934416525e-04, -8.91196859397070326e-06};
    const double z = 0.25 * x * x;
    const double i1 =
        (((rational(p_i, q_i, z) + 0.08695471286773681640625) * z * z +
          0.5 * z + 1.0) *
         x) *
        0.5;

    // K_1(x) = 1/x + x*R(x*x) + log(x)*I_1(x).
    static const double p_k[] = {
        -3.07965757829206184e-01, -7.80929703673074907e-02,
        -2.70619343754051620e-03, -2.49549522229072008e-05};
    static const double q_k[] = {
        1.00000000000000000e+00, -2.36316836412163098e-02,
        2.64524577525962719e-04, -1.49749618004162787e-06};
    const double k1 =
        rational(p_k, q_k, x * x) * x + 1.0 / x + std::log(x) * i1;
    // Safe here: 0 < x <= 1, so exp(x) cannot overflow.
    return std::exp(x) * k1;
  }

  // For x > 1, approximate sqrt(x)*exp(x)*K_1(x) directly as
  // Y + P(1/x)/Q(1/x). No exponential is evaluated in this branch.
  static const double p[] = {-1.97028041029226295e-01, -2.32408961548087617e+00,
                             -7.98269784507699938e+00, -2.39968410774221632e+00,
                             3.28314043780858713e+01,  5.67713761158496058e+01,
                             3.30907788466509823e+01,  6.62582288933739787e+00,
                             3.08851840645286691e-01};
  static const double q[] = {1.00000000000000000e+00, 1.41811409298826118e+01,
                             7.35979466317556420e+01, 1.77821793937080859e+02,
                             2.11014501598705982e+02, 1.19425262951064454e+02,
                             2.88448064302447607e+01, 2.27912927104139732e+00,
                             2.50358186953478678e-02};
  return (1.4503421783447265625 + rational(p, q, 1.0 / x)) / std::sqrt(x);
}

/*
Copyright (c) 2006 Xiaogang Zhang
Copyright (c) 2017 John Maddock

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
