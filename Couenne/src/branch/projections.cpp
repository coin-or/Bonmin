/*
 * Name:    projections.c
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tools for projecting points on lines/planes
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <stdio.h>
#include <stdlib.h>

#include <CouenneTypes.hpp>
#include <CouennePrecisions.hpp>

/*  compute projection of point (x0, y0) on the segment defined by
 *  line ax + by + c <>= 0 (sign provided by parameter sign) and
 *  bounds [lb, ub] on x. Return distance from segment, 0 if satisfied
 */

CouNumber project (CouNumber a, CouNumber b, CouNumber c, 
		   CouNumber x0, CouNumber y0, 
		   CouNumber lb, CouNumber ub, int sign,
		   CouNumber *xp, CouNumber *yp) {

  /* compute projection of (x0,y0) onto line ax+by+c=0 */
  register CouNumber
    t  = - (a*x0 + b*y0 + c);

  /* projection coordinates */
  CouNumber xpr, ypr;

  /* does point lie on line? */
  if (fabs (t) < COUENNE_EPS) return 0.; 

  /* check if point satisfies inequality */
  if      (sign > 0) {if (t < 0.) return 0.;}
  else if (sign < 0) {if (t > 0.) return 0.;}

  /* t corresponding to intersection point */
  t /= sqrt (a*a + b*b);

  /* compute projection coordinates */
  xpr = x0 + a*t;
  ypr = y0 + b*t;

  /* don't need sign any longer, take its absolute value */
  if (t < 0.) t = -t;

  /* if projected point is outside [lb,ub], set xp to closest bound
     and yp accordingly, and compute distance to (x0,y0) */
  if ((xpr < lb) || (xpr > ub)) {

    if      (xpr < lb) xpr = lb; 
    else if (xpr > ub) xpr = ub; 

    ypr = (- c - a * xpr) / b - y0;
    xpr -= x0;

    t = sqrt (xpr * xpr + ypr * ypr);
  }

  /* update output parameters */
  if (xp) *xp = xpr;
  if (yp) *yp = ypr;

  /* return distance */
  return t;
}


/*  Compute best branching point (within interval [lb,ub]) given
 *  function f, its derivative fp and the current optimum (x0,y0)
 */
/*
CouNumber bestBranchPoint (CouNumber lb, CouNumber ub, 
			   CouNumber x0, CouNumber y0, 
			   unary_function f,
			   unary_function fp) {

  return 0;
}
*/
/*
int main (int argc, char **argv) {

  CouNumber 
    a = atof (argv [1]),
    b = atof (argv [2]),
    c = atof (argv [3]),
    x0 = atof (argv [4]),
    y0 = atof (argv [5]),
    lb = atof (argv [6]),
    ub = atof (argv [7]);

  int sign = atoi (argv [8]);

  char sig = (sign < 0) ? '<' : (sign > 0) ? '>' : ' ';

  printf ("projecting (%.3f,%.3f) on %.3f x + %.3f y + %.3f %c= 0, xp in [%.3f,%.3f]\n",
	  x0, y0, a, b, c, sig, lb, ub);

  printf (" ==> distance is %.4f\n",
	  project (a, b, c, x0, y0, lb, ub, sign));

}
*/
