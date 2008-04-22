/*
 * Name:    branchExprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for exponentials
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprExp.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "projections.hpp"
#include "funtriplets.hpp"


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprExp::selectBranch (const CouenneObject *obj, 
				 const OsiBranchingInformation *info,
				 expression *&var, 
				 double * &brpts, 
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {

  // two cases: inside and outside the curve. 
  //
  // Inside: the distance depends on the projection of the current
  // point onto the would-be upper envelopes, which forces us to look
  // at it numerically. If both bounds are infinite, create a ThreeWay
  // branch.
  //
  // Outside: it suffices to project the current point on the line
  // (i.e. to get the closest point on the line) to get the maxi-min
  // displacement.
  //
  // As for all monotone functions, after choosing *brpts it is
  // equivalent to choose w's or x's index as ind, as the implied- and
  // propagated bounds will do the rest.

  var = argument_;

  brDist = (double *) realloc (brDist, 2 * sizeof (double));
  brpts  = (double *) realloc (brpts, sizeof (double));

  int
    ind = var -> Index (),
    wi  = obj -> Reference () -> Index ();

  assert ((ind >= 0) && (wi >= 0));

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  // Outside //////////////////////////////////////////////////////////////////

  if (y0 < exp (x0)) {

    // Look for point (x1,y1) on curve y=exp(x) closest to current
    // (x0,y0), branch at x1

    *brpts = obj -> midInterval (powNewton (x0, y0, exp, exp, exp), l, u);

    way = TWO_RAND;

    y0 -= exp (*brpts);
    x0 -= *brpts;

    return sqrt (brDist [0] = brDist [1] = sqrt (x0*x0 + y0*y0)); // exact distance
  }

  // Inside. Four cases: ///////////////////////////////////////////////////////

  if ((l < -COUENNE_INFINITY) && 
      (u >  COUENNE_INFINITY)) { // unbounded in both directions

    /*    // TODO: restore when we can do three-way branching
#if 0
    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first
    brpts [0] = x0;       // draw vertical   from (x0,y0) south to curve y=exp(x)
    brpts [1] = log (y0); //      horizontal              east
    CouNumber 
      a = y0 - exp (x0),  // sides of a triangle with (x0,y0)
      b = log (y0) - x0;  // as one of the vertices
    // exact distance from central interval, from others it's a and b
    return (a * cos (atan (a/b))); 
#endif    */

    // follow South-East diagonal to find point on curve
    // so that current point is surely cut 
    *brpts = 0.5 * (x0 + log (y0)); 
    way = TWO_RAND;

    return CoinMin (brDist [0] = log (y0) - x0, 
		    brDist [1] = y0 - exp (x0));
  }

  // 2,3) at least one of them is finite

  if (l < - COUENNE_INFINITY) { // 2) unbounded from below --> break vertically

    *brpts = obj -> midInterval (x0, l, u);

    way = TWO_RIGHT;
    return CoinMin (brDist [0] = y0 - exp (x0), 
		    brDist [1] = projectSeg (x0, y0, *brpts, exp (*brpts), u, exp (u), -1));
  } 

  if (u > COUENNE_INFINITY) { // 3) unbounded from above -- break horizontally

    *brpts = obj -> midInterval (log (y0), l, u);

    way = TWO_LEFT;
    return CoinMin (brDist [0] = log (y0) - x0, 
		    brDist [1] = projectSeg (x0, y0, l, exp (l), *brpts, exp (*brpts), -1));
  }

  // 4) both are finite

  simpletriplet ft (exp, exp, exp, log);
  *brpts = obj -> getBrPoint (&ft, x0, l, u); // select based on strategy

  way = TWO_RAND;

  // exact distance
  return CoinMin (brDist [0] = projectSeg (x0, y0, l, exp (l), *brpts, exp (*brpts),             -1),
		  brDist [1] = projectSeg (x0, y0,             *brpts, exp (*brpts), u, exp (u), -1));
}
