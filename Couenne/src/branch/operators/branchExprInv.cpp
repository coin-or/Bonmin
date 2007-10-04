/*
 * Name:    branchExprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch selection for 1/x
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprInv.hpp>
#include <exprDiv.hpp>
#include <exprPow.hpp>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments

CouNumber exprInv::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {

  // two cases: inside or outside the bellies. 
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

  // two cases: 
  // 
  // 1) bounds include 0: three way branch to exclude 0 (refer to branchExprDiv.cpp)
  // 2) otherwise         two   way branch

  if ((l < 0) && (u > 0)) {

    brpts = (double *) realloc (brpts, 2 * sizeof (CouNumber));
    way = THREE_RIGHT;

    brpts [0] = (l >= - BR_NEXT_ZERO - COUENNE_EPS) ? (l * BR_MULT) : - BR_NEXT_ZERO;
    brpts [1] = (u <=   BR_NEXT_ZERO + COUENNE_EPS) ? (u * BR_MULT) :   BR_NEXT_ZERO;

    return ((fabs (x0) < COUENNE_EPS) ? 1. : 
	    fabs (1 / x0 - info -> solution_ [wi]));
  }

  // case 2: look if inside or outside of belly (refer to branchExprExp.cpp)

  if (x0*y0 < 1) { // outside bellies 

    brpts = (double *) realloc (brpts, sizeof (double));

    *brpts = powNewton (x0, y0, inv, oppInvSqr, inv_dblprime);

    if      (*brpts < l) *brpts = (u <   COUENNE_INFINITY) ? ((l+u)/2) : l + 1;
    else if (*brpts > u) *brpts = (l > - COUENNE_INFINITY) ? ((l+u)/2) : u - 1;

    way = (u < 0) ? TWO_RIGHT : TWO_LEFT; // explore finite interval first

    CouNumber dy = y0 - 1. / *brpts;
    x0 -= *brpts;
    return sqrt (x0*x0 + dy*dy);
  }

  // Inside, x0*y0 >= 1. Two cases:
 
  if ((l <   COUENNE_EPS) && (u >   COUENNE_INFINITY) || 
      (u > - COUENNE_EPS) && (l < - COUENNE_INFINITY)) {

    // 1) bounds are infinite both horizontally and vertically
    //    (i.e. [-inf,0] or [0,inf]) --> three way branching

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = x0;      // draw vertical   from (x0,y0) south (north) to curve y=1/x
    brpts [1] = 1. / y0; //      horizontal              west  (east)

    CouNumber a = fabs (y0 - 1 / x0), // sides of a triangle with (x0,y0)
              b = fabs (x0 - 1 / y0), // as one of the vertices
              c = a * cos (atan (a/b));

    return mymin (a, mymin (b, c));
  }

  // 2) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l < - COUENNE_INFINITY) { // and u away from -0

    *brpts = x0;
    if (*brpts > u - COUENNE_NEAR_BOUND) *brpts = u-1;
    way = TWO_RIGHT;
    return y0 - 1. / x0;
  } 

  if (u > COUENNE_INFINITY) { // and l away from +0

    *brpts = x0;
    if (*brpts < l + COUENNE_NEAR_BOUND) *brpts = l+1;
    way = TWO_LEFT;
    return y0 - 1. / x0;
  }

  // last case: nice finite interval and limited curve

  *brpts = powNewton (x0, y0, inv, oppInvSqr, inv_dblprime);

  if ((*brpts > u - COUENNE_NEAR_BOUND) ||
      (*brpts < l + COUENNE_NEAR_BOUND))
    *brpts = (l+u) / 2;

  way = TWO_RAND;

  CouNumber dx = x0 - *brpts,
    dy = y0 - 1. / *brpts;

  return sqrt (dx*dx + dy*dy);
}
