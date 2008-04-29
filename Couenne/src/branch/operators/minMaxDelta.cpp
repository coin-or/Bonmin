/*
 * Name:    minMaxDelta.cpp
 * Author:  Pietro Belotti
 * Purpose: general function for computing best branching point based
 *          on min max height of resulting convexifications (dychotomic 
 *          search)
 *
 * (C) Carnegie-Mellon University, 2007-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneObject.hpp"
#include "funtriplets.hpp"

const int maxIter = 20;

///
CouNumber curvDistance (funtriplet *ft, CouNumber lb, CouNumber ub) {

  // Consider the function f(x) between lb and ub. The slope of the
  // convexification on the concave side, y = alpha x + alpha0, is:

  CouNumber alpha = (ft -> F (ub) - ft -> F (lb)) / (ub - lb);

  // and the constant term, alpha0, is

  CouNumber alpha0 = (ub * ft -> F (lb) - lb * ft -> F (ub)) / (ub - lb);

  // The point at which f(.) has derivative equal to the slope is the
  // point of maximum height w.r.t the slope. The point z where
  // maximum of f(z) - (ax+b), where (ax+b) is the convexification
  // line (a=slope), is such that 
  //
  // f'(z) - alpha = 0   ==> z = (f')^{-1} (alpha)

  CouNumber z = ft -> FpInv (alpha);

  // The real height is computed as [f(z) - (alpha * z + alpha0)]
  // divided by the norm sqrt (alpha^2 + 1)

  return ((ft -> F (z) - (alpha * z + alpha0)) / sqrt (alpha * alpha + 1));
}


///
CouNumber minMaxDelta (funtriplet *ft, CouNumber lb, CouNumber ub) {

  CouNumber 
    lbm = lb,                // extremes of the interval where to look 
    ubm = ub,     
    b   = 0.5 * (lbm + ubm); // starting point

  for (int iter = 0; iter < maxIter; iter++) {

    CouNumber distL = curvDistance (ft, lb,  b),  // max height at left
              distR = curvDistance (ft,  b, ub),  // max height at right
              delta = fabs (distL) - fabs (distR);

    //    fprintf (stderr, "%4d %10g %10g %10g %10g %10g %10g\n", 
    //	     iter, lbm, ubm, b, distL, distR, delta);

    if (fabs (delta) < COUENNE_EPS) 
      break;

    CouNumber oldb = b;

    // choose a smarter b based on an estimate of the derivative of
    // the distance function at the current point, knowing it's null
    // at the extremes

    b = 0.5 * (lbm + ubm);

    if (delta > 0) ubm = oldb; // right max height is smaller, move left
    else           lbm = oldb; // and viceversa
  }

  return b;
  //return obj -> midInterval (b, lb, ub);
}


///
CouNumber maxHeight (funtriplet *ft, CouNumber lb, CouNumber ub) {
  /* fprintf (stderr,"slope is (%g - %g) / (%g - %g) = %g / %g = %g ----> inverse is %g\n", 
	  ft -> F (ub), 
	  ft -> F (lb), 
	  ub, lb,
	  ft -> F (ub) - ft -> F (lb),
	  (ub - lb),
	  (ft -> F (ub) - ft -> F (lb)) / (ub - lb),
	  ft -> FpInv ((ft -> F (ub) - ft -> F (lb)) / (ub - lb)));*/
  return (ft -> FpInv ((ft -> F (ub) - ft -> F (lb)) / (ub - lb)));
}
