/*
 * Name:    projections.h
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tools for projecting points on lines/planes
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouennePrecisions.h>

extern "C" {

  /*  compute projection of point (x0, y0) on the segment defined by
   *  line ax + by + c <>= 0 (sign provided by parameter sign) and
   *  bounds [lb, ub] on x. Return distance from segment, 0 if satisfied
   */

  CouNumber project (CouNumber a, CouNumber b, CouNumber c, 
		     CouNumber x0, CouNumber y0, 
		     CouNumber lb, CouNumber ub, int sign,
		     CouNumber *, CouNumber *);


  /*  Compute best branching point (within interval [lb,ub]) given
   *  function f, its derivative fp and the current optimum (x0,y0)
   */
  CouNumber bestBranchPoint (CouNumber lb, CouNumber ub, 
			     CouNumber x0, CouNumber y0, 
			     unary_function f,
			     unary_function fp);
}
