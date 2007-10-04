/*
 * Name:    trigNewton.cpp
 * Author:  Pietro Belotti
 * Purpose: numerically find tangents to (co)sines
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <CouenneTypes.hpp>

#define MAX_ITER 1000
#define COU_TRIG_TOLERANCE 1e-12

CouNumber trigNewton (CouNumber a, CouNumber l, CouNumber u) {

  // Find a zero to the function
  //
  // F(x) = cos x - (sin x - sin a) / (x - a)
  // 
  // whose derivative is
  //
  // F'(x) = - sin x - cos x / (x-a) + (sin x - sin a) / (x - a)^2

  if (l>u) {
    register CouNumber swap = l;
    l = u;
    u = swap;
  }

  register CouNumber xk = 0.5 * (u+l);

  CouNumber sina  = sin (a),
            sinxk = sin (xk),
            cosxk = cos (xk),
            dy    = sinxk - sina,
            dx    = xk - a,
            dydx  = dy/dx,
            F     = cosxk - dydx;

  // Newton loop. Tolerance is set above
  for (register int k = MAX_ITER; (fabs (F) > COU_TRIG_TOLERANCE) && k--;) {

    CouNumber Fp = sinxk + (cosxk - dydx) / dx;

    xk += F/Fp;

    if      (xk < l) xk = l;
    else if (xk > u) xk = u;

    sinxk = sin (xk);
    cosxk = cos (xk);
    dy    = sinxk - sina;
    dx    = xk - a;
    dydx  = dy/dx;
    F     = cosxk - dydx;
  }

  return xk;
}

/*
int main (int argc, char **argv) {

  CouNumber a = atof (argv [1]),
            l = atof (argv [2]),
            u = atof (argv [3]), r;

  for (register int i=100000; i--;)
    r = trigNewton (a, l, u);

  printf ("b0 = %.14f: slope %.15f, derivative %.15f\n", 
	  r, (sin (r) - sin (a)) / (r-a), cos (r));

  return 0;
}
*/
