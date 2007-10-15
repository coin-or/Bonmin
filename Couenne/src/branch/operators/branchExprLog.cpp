/*
 * Name:    branchExprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for logarithms
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprLog.hpp"
#include "exprPow.hpp"

#include "CoinHelperFunctions.hpp"

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "projections.hpp"


#define SQ_COUENNE_EPS COUENNE_EPS * COUENNE_EPS

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprLog::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {

  // quite similar to exprExp::selectBranch() (see branchExprExp.cpp)

  // two cases: inside or outside the belly. 
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
  // As for all monotonous functions, after choosing *brpts it is
  // equivalent to choose w's or x's index as ind, as the implied- and
  // propagated bounds will do the rest.

  ind    = argument_ -> Index ();
  int wi = w         -> Index ();

  assert ((ind >= 0) && (wi >= 0));

  CouNumber 
    y0 = info -> solution_ [wi],
    x0 = info -> solution_ [ind],
    l  = info -> lower_    [ind],
    u  = info -> upper_    [ind];

  if (u < COUENNE_EPS) { // strange case, return default branching rule
    ind = -1; 
    return 0;
  }

  if (x0 < SQ_COUENNE_EPS) // very unlikely...
    x0 = SQ_COUENNE_EPS;

  if (y0 > log (x0)) { 

    // Outside

    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = midInterval (powNewton (x0, y0, log, inv, oppInvSqr), l, u);

    /*    if      (*brpts < l) if (u <   COUENNE_INFINITY) *brpts = (l+u)/2;
                         else                        *brpts = 10*l + 1;
    else if (*brpts > u) if (l < - COUENNE_INFINITY) *brpts = u/2;
                         else                        *brpts = (l+u)/2; */

    way = TWO_LEFT;
    CouNumber dy = y0 - log (*brpts);
    x0 -= *brpts;

    return sqrt (x0*x0 + dy*dy); // exact distance
  } 

  // Inside. Two cases:
 
  if ((l <= SQ_COUENNE_EPS) && 
      (u > COUENNE_INFINITY)) {

    // 1) curve is unlimited in both senses --> three way branching

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = exp (y0); // draw horizontal from (x0,y0) east  to curve y=log(x)
    brpts [1] = x0;       //      vertical                north

    CouNumber a = x0 - exp (y0), // sides of a triangle with (x0,y0)
              b = log (x0) - y0, // as one of the vertices
              c = a * cos (atan (a/b));

    return CoinMin (a, CoinMin (b, c)); // exact distance
  } 

  // 2) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l <= SQ_COUENNE_EPS) { // u is finite

    *brpts = midInterval (exp (y0), l, u);

    /*if ((*brpts > u - COUENNE_NEAR_BOUND) ||
	(*brpts < l + COUENNE_NEAR_BOUND)) 
	*brpts = (l+u) / 2;*/

    way = TWO_RIGHT;
    return projectSeg (x0, y0, *brpts, log (*brpts), x0, log (x0), +1); // exact distance

    //    return CoinMin (x0 - exp (y0), log (x0) - y0);
  }
 
  if (u > COUENNE_INFINITY) { // l is far from zero

    *brpts = midInterval (x0, l, u);
    way = TWO_LEFT;

    /*if (*brpts < l + COUENNE_NEAR_BOUND)
     *brpts = l+1;*/

    return projectSeg (x0, y0, *brpts, log (*brpts), x0, log (x0), +1); // exact distance
    //return log (x0) - y0;
  } 

  // both are finite

  //  *brpts = midInterval (powNewton (x0, y0, log, inv, oppInvSqr), l, u); 
  // WRONG! Local minima may be at bounds

  *brpts = midInterval (x0, l, u); 

  /*if ((*brpts > u - COUENNE_NEAR_BOUND) ||
      (*brpts < l + COUENNE_NEAR_BOUND))
      *brpts = (l+u) / 2;*/

  way = TWO_RAND;

  CouNumber 
    dx = x0 - *brpts,
    dy = y0 - log (*brpts);

  return sqrt (dx*dx + dy*dy); // exact distance
}
