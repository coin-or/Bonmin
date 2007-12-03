/*
 * Name:    powNewton.cpp
 * Author:  Pietro Belotti
 * Purpose: numerically find tangents to power functions
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "CouenneTypes.hpp"

#define MAX_ITER 10
#define COU_POW_TOLERANCE 1e-12

//#define DEBUG_POWNEW

#ifndef DEBUG_POWNEW
#include "funtriplets.hpp"
#else
#include <stdlib.h>
#include <stdio.h>
#endif

CouNumber powNewton (CouNumber xc, CouNumber yc, 
		     unary_function f, 
		     unary_function fp, 
		     unary_function fpp) {

  // Find a zero to the function
  //
  // F(x) = x - xc + f'(x) (f(x) - yc)
  //
  // where f(x) is either x^k, exp(x), or log(x).
  // The derivative of F(x) is
  //
  // F'(x) = 1 + f''(x) (f(x) - yc) + (f'(x))^2
  //
  // Apply usual update:
  //
  // x(k+1) = x(k) - F(x(k))/F'(x(k))

  register CouNumber xk = xc;

  CouNumber fk  = f (xk) - yc,
            fpk = fp (xk),
            F   = fpk * fk,
            Fp  = 1 + fpp (xk) * fk + fpk * fpk;

  // Newton loop. Tolerance is set above
  for (int k = MAX_ITER; k--;) {

    xk -= F / Fp;

    fk  = f (xk) - yc;
    fpk = fp (xk);
    F   = xk - xc + fpk * fk;

    //    printf ("xk = %g; F = %g, fk = %g, fpk = %g\n", xk, F, fk, fpk);

    if (fabs (F) < COU_POW_TOLERANCE) break;
    Fp  = 1 + fpp (xk) * fk + fpk * fpk;
  }

  return xk;
}

#ifndef DEBUG_POWNEW

///
CouNumber powNewton (CouNumber xc, CouNumber yc, funtriplet *tri) {

  // Find a zero to the function
  //
  // F(x) = x - xc + f'(x) (f(x) - yc)
  //
  // where f(x) is either x^k, exp(x), or log(x).
  // The derivative of F(x) is
  //
  // F'(x) = 1 + f''(x) (f(x) - yc) + (f'(x))^2
  //
  // Apply usual update:
  //
  // x(k+1) = x(k) - f(x(k))/f'(x(k))

  register CouNumber xk = xc;

  CouNumber fk  = tri -> F (xk) - yc,
            fpk = tri -> Fp (xk),
            F   = fpk * fk,
            Fp  = 1 + tri -> Fpp (xk) * fk + fpk * fpk;

  // Newton loop. Tolerance is set above
  for (int k = MAX_ITER; k--;) {

    xk -= F / Fp;

    fk  = tri -> F (xk) - yc;
    fpk = tri -> Fp (xk);
    F   = xk - xc + fpk * fk;
    if (fabs (F) < COU_POW_TOLERANCE) break;
    Fp  = 1 + tri -> Fpp (xk) * fk + fpk * fpk;
  }

  return xk;
}
#else

/// the operator itself
inline CouNumber inv (register CouNumber arg) 
{return 1.0 / arg;}


/// derivative of inv (x)
inline CouNumber oppInvSqr (register CouNumber x) 
{return (- inv (x*x));}


/// inv_dblprime, second derivative of inv (x)
inline CouNumber inv_dblprime (register CouNumber x) 
{return (2 * inv (x*x*x));}


int main (int argc, char **argv) {

  CouNumber r, 
    xc = atof (argv [2]),
    yc = atof (argv [3]);

  unary_function 
    f   = log,
    fp  = inv,
    fpp = oppInvSqr;

  //expon = atof (argv [1]);

  for (register int i=1; i--;)
    r = powNewton (xc, yc, f, fp, fpp);

  printf ("xc = %.14f: xk = %.15f, slope %.15f -- %.15f ==> [%.15f = -1?]\n", 
	  xc, r, fp (r), 
	           (yc - f (r)) / (xc - r), 
	  fp (r) * (yc - f (r)) / (xc - r));

  return 0;
}

#endif
