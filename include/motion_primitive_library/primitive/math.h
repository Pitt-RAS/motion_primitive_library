/**
 * @file math.h
 * @brief Polynomial roots solver

 * Solving real roots for n-th order polynomial:
    if n < 5, the closed form solution will be calculated;
    if n >= 5, using Eigen Polynomials solver which is slower but correct.
 */
#pragma once
#include <motion_primitive_library/common/data_type.h>

// BAD HEADERS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <unsupported/Eigen/Polynomials>
#pragma GCC diagnostic pop
// END BAD HEADERS

#include <iostream>

///Cubic equation: \f$a*t^3+b*t^2+c*t+d = 0\f$
int cubic(decimal_t a,
          decimal_t b,
          decimal_t c,
          decimal_t d,
          std::array<decimal_t, 4>& out);

///Quartic equation: \f$a*t^4+b*t^3+c*t^2+d*t+e = 0\f$
int quartic(decimal_t a,
            decimal_t b,
            decimal_t c,
            decimal_t d,
            decimal_t e,
            std::array<decimal_t, 4>& out);

/*! \brief General solver for \f$a*t^4+b*t^3+c*t^2+d*t+e = 0\f$

  \f$a, b, c\f$ can be zero. The function itself checks the highest order of the polynomial.
  */
std::vector<decimal_t> solve(decimal_t a, decimal_t b, decimal_t c, decimal_t d, decimal_t e);

int solve_fast(decimal_t a,
               decimal_t b,
               decimal_t c,
               decimal_t d,
               decimal_t e,
               std::array<decimal_t, 4>& out);

///A more general solver for \f$a*t^6+b*t^5+c*t^4+d*t^3+e*t^2+f*t+g = 0\f$
std::vector<decimal_t> solve(decimal_t a, decimal_t b, decimal_t c, decimal_t d, decimal_t e, decimal_t f, decimal_t g);
int solve_fast(decimal_t a,
               decimal_t b,
               decimal_t c,
               decimal_t d,
               decimal_t e,
               decimal_t f,
               decimal_t g,
               std::array<decimal_t, 6>& out);

///Return \f$n!\f$
int factorial(int n);

///Return \f$t^n\f$
decimal_t power(decimal_t t, int n);
 
