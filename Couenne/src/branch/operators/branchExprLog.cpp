/*
 * Name:    branchExprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for logarithms
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprLog.hpp>
#include <exprPow.hpp>

#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


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

  if ((ind < 0) || (wi < 0)) {printf ("Couenne, w=log(x): negative index\n"); exit (-1);}

  CouNumber y0 = info -> solution_ [wi],
    x0 = info -> solution_ [ind],
    l  = info -> lower_    [ind],
    u  = info -> upper_    [ind];

  if (y0 > log (x0)) { 

    if (x0 == 0) x0 = 1e-30;

    // Outside

    brpts = (double *) realloc (brpts, sizeof (double));

    *brpts = powNewton (x0, y0, log, inv, oppInvSqr);

    if      (*brpts < l) *brpts = (u < COUENNE_INFINITY) ? ((l+u)/2) : (10*l + 1);
    else if (*brpts > u) *brpts = (l+u) / 2;

    way = TWO_LEFT;
    CouNumber dy = y0 - log (*brpts);
    x0 -= *brpts;
    return sqrt (x0*x0 + dy*dy);

  } 

  // Inside. Two cases:
 
  if ((l < COUENNE_EPS * COUENNE_EPS) && 
      (u > COUENNE_INFINITY)) {

    // 1) curve is unlimited in both senses --> three way branching

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = exp (y0); // draw horizontal from (x0,y0) south to curve y=log(x)
    brpts [1] = x0;       //      vertical                east

    CouNumber a = x0 - exp (y0), // sides of a triangle with (x0,y0)
      b = y0 - log (x0), // as one of the vertices
      c = a * cos (atan (a/b));

    return mymin (a, mymin (b, c));

  } 

  // 2) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l < COUENNE_EPS * COUENNE_EPS) {

    *brpts = exp (y0);
    if ((*brpts > u - COUENNE_NEAR_BOUND) ||
	(*brpts < l + COUENNE_NEAR_BOUND)) 
      *brpts = (l+u) / 2;

    way = TWO_RIGHT;
    return mymin (x0 - exp (y0), y0 - log (x0));

  }
 
  if (u > COUENNE_INFINITY) {

    *brpts = x0;
    if (*brpts < l + COUENNE_NEAR_BOUND) *brpts = l+1;
    way = TWO_LEFT;
    return y0 - log (x0);

  } 

  // both are finite

  *brpts = powNewton (x0, y0, log, inv, oppInvSqr);

  if ((*brpts > u - COUENNE_NEAR_BOUND) ||
      (*brpts < l + COUENNE_NEAR_BOUND))
    *brpts = (l+u) / 2;

  way = TWO_RAND;

  CouNumber 
    dx = x0 - *brpts,
    dy = y0 - log (*brpts);

  return sqrt (dx*dx + dy*dy);
}
