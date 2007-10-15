/*
 * Name:    branchExprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for exponentials
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprExp.hpp"
#include "exprPow.hpp"

#include "CoinHelperFunctions.hpp"

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "projections.hpp"


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprExp::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {

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

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  if (y0 < exp (x0)) {

    // Outside: look for point (x1,y1) on curve y=exp(x) closest to
    // current (x0,y0), branch at x1

    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = midInterval (powNewton (x0, y0, exp, exp, exp), l, u);

    //    if      (*brpts < l) *brpts = (u <   COUENNE_INFINITY) ? ((l+u)/2) : +1 + (l>0) ? l*2 : l/2;
    //    else if (*brpts > u) *brpts = (l > - COUENNE_INFINITY) ? ((l+u)/2) : -1 + (u<0) ? u*2 : u/2;

    way = TWO_RAND;
    CouNumber dy = y0 - exp (*brpts);
    x0 -= *brpts;

    return sqrt (x0*x0 + dy*dy); // exact distance
  }

  // Inside. Four cases: ///////////////////////////////////////////////////////

  // 1) both bounds are infinite --> three way branching
 
  if ((l < -COUENNE_INFINITY) && 
      (u >  COUENNE_INFINITY)) {

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = x0;       // draw vertical   from (x0,y0) south to curve y=exp(x)
    brpts [1] = log (y0); //      horizontal              east

    CouNumber 
      a = y0 - exp (x0),        // sides of a triangle with (x0,y0)
      b = x0 - log (y0),        // as one of the vertices
      c = a * cos (atan (a/b)); // all three quantities are nonnegative

    return CoinMin (a, CoinMin (b, c)); // exact distance
  }

  // 2,3,4) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l < -COUENNE_INFINITY) { // 2) unbounded from below

    *brpts = midInterval (x0, l, u);
    /*    if (*brpts > u - COUENNE_NEAR_BOUND) 
     *brpts = u-1;*/

    way = TWO_RIGHT;
    return CoinMin (y0 - exp (x0), projectSeg (x0, y0, x0, exp (x0), u, exp (u), -1));
  } 

  if (u > COUENNE_INFINITY) { // 3) unbounded from above

    *brpts = midInterval (log (y0), l, u);

    /*    if (*brpts < l + COUENNE_NEAR_BOUND) 
     *brpts = l+1;*/

    way = TWO_LEFT;
    return CoinMin (x0 - log (y0), projectSeg (x0, y0, l, exp (l), log (y0), y0, -1));
  }

  // 4) both are finite

  // find closest point on curve
  *brpts = midInterval (powNewton (x0, y0, exp, exp, exp), l, u);

  /*
    if ((*brpts > u - COUENNE_NEAR_BOUND) ||
        (*brpts < l + COUENNE_NEAR_BOUND))
    *brpts = (l+u) / 2;
    */

  way = TWO_RAND;

  CouNumber dx = x0 - *brpts,
    dy = y0 - exp (*brpts);

  return sqrt (dx*dx + dy*dy);
}
