/*
 * Name:    powNewton.cpp
 * Author:  Pietro Belotti
 * Purpose: numerically find tangents to power functions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>
#include <CouenneTypes.h>

#define MAX_ITER 1000
#define COU_POW_TOLERANCE 1e-12

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
  // Follow usual update:
  //
  // x(k+1) = x(k) - f(x(k))/f'(x(k))

  register CouNumber xk = xc;

  CouNumber fk  = f (xk) - yc,
            fpk = fp (xk),
            F   = fpk * fk,
            Fp  = 1 + fpp (xk) * fk + fpk * fpk;

  // Newton loop. Tolerance is set above
  for (register int k = MAX_ITER; (fabs (F) > COU_POW_TOLERANCE) && k--;) {

    xk -= F / Fp;

    fk  = f (xk) - yc;
    fpk = fp (xk);
    F   = xk - xc + fpk * fk;
    Fp  = 1 + fpp (xk) * fk + fpk * fpk;
  }

  return xk;
}

/*
CouNumber expon = 4;

inline CouNumber f (CouNumber x)
{return pow (x, expon);}

inline CouNumber fp (CouNumber x) 
{return expon * pow (x, expon-1);}

inline CouNumber fpp (CouNumber x) 
{return expon * (expon-1) * pow (x, expon-2);}

int main (int argc, char **argv) {

  CouNumber r, 
    xc = atof (argv [2]),
    yc = atof (argv [3]);

  expon = atof (argv [1]);

  for (register int i=10000; i--;)
    r = powNewton (xc, yc, f, fp, fpp);

  printf ("xc = %.14f: xk = %.15f, slope %.15f -- %.15f ==> %.15f\n", 
	  xc, r, fp (r), (yc - f (r)) / (xc - r), fp (r) * (yc - f (r)) / (xc - r));

  return 0;
}
*/
