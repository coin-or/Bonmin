/*
 * Name:    branchExprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: select branch for divisions
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprDiv.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneObject.hpp"


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprDiv::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index (),
      wi = w            -> Index ();

  assert ((xi >= 0) && (yi >= 0));

  // choosing branching variable and -point is difficult, use
  // proportion in bound intervals

  CouNumber yl = info -> lower_    [yi], 
            yu = info -> upper_    [yi],
            y0 = info -> solution_ [yi];

  // if [yl,yu] contains 0, create three nodes, with a small central
  // one tight around zero -- it has a bad convexification and will be
  // explored last

  if ((yl < 0) && (yu > 0)) {

    ind = yi;
    //    brpts = (double *) realloc (brpts, 2 * sizeof (CouNumber));
    //    way = THREE_RIGHT;
    way = TWO_RAND;
    brpts = (double *) realloc (brpts, sizeof (CouNumber));

    brpts [0] = 0;
    //    brpts [0] = (yl >= -BR_NEXT_ZERO - COUENNE_EPS) ? (yl * BR_MULT) : -BR_NEXT_ZERO;
    //    brpts [1] = (yu <=  BR_NEXT_ZERO + COUENNE_EPS) ? (yu * BR_MULT) :  BR_NEXT_ZERO;

    return ((fabs (y0) < COUENNE_EPS) ? 1. : 
	    fabs (info -> solution_ [xi] / y0 - info -> solution_ [w -> Index ()]));
  }

  // From now on, [yl,yu] may be unlimited in one sense only, and
  // interval does not contain 0.
  //
  // As convexification is still carried out by applying McCormick
  // rules to x=w*y (where original operator is w=x/y), try to get
  // closer to a point where both y and w are bounded, if necessary by
  // branching on w.
  //
  // First, branch on y if unbounded, and then on w. As a result of
  // bound tightening, if both y and w are bounded, x is, too.

  if ((yl < -COUENNE_INFINITY) ||
      (yu >  COUENNE_INFINITY)) {

    ind = yi;
    brpts = (double *) realloc (brpts, sizeof (CouNumber));

    // if y0 close to bounds, branch away from it
    if      (fabs (y0-yl) < COUENNE_NEAR_BOUND) *brpts = y0 + 1 + yl*10;
    else if (fabs (y0-yu) < COUENNE_NEAR_BOUND) *brpts = y0 - 1 + yu*10;
    else                                        *brpts = y0;

    way = (y0 > 0) ? TWO_LEFT : TWO_RIGHT;

    return ((fabs (y0) < COUENNE_EPS) ? 1. : 
	    fabs (info -> solution_ [xi] / y0 - info -> solution_ [w -> Index ()]));
  }

  // y is bounded, and y0 should not be 0; if w is unbounded, it is
  // better to bound w first as x will be too.

  CouNumber wl = info -> lower_    [wi], 
            wu = info -> upper_    [wi],
            w0 = info -> solution_ [wi],
            x0 = info -> solution_ [xi];

  if ((wl < -COUENNE_INFINITY) || 
      (wu >  COUENNE_INFINITY)) {

    ind = wi;

    if ((wl < -COUENNE_INFINITY) &&
	(wu >  COUENNE_INFINITY)) {

      // unbounded in two directions

      CouNumber wreal = x0 / y0, 
	wmin = w0, wmax = wreal; // (x0,y0,w0) is below surface

      if (wreal < w0) { // (x0,y0,w0) is above surface
	wmin = wreal;
	wmax = w0;
      }

      // unbounded in both directions: three way branching 
      brpts = (double *) realloc (brpts, 2 * sizeof (CouNumber));
      brpts [0] = wmin;
      brpts [1] = wmax;
      way = THREE_CENTER;

    } else {

      // unbounded in one direction only, use two way branching

      brpts = (double *) realloc (brpts, sizeof (CouNumber));

      // if y0 close to bounds, branch away from it
      if      (fabs (w0-wl) < COUENNE_NEAR_BOUND) *brpts = w0 + 1 + wl*10;
      else if (fabs (w0-wu) < COUENNE_NEAR_BOUND) *brpts = w0 - 1 + wu*10;
      else                                        *brpts = w0;

      way = (wl < - COUENNE_INFINITY) ? TWO_RIGHT : TWO_LEFT;
    }

    return ((fabs (y0) < COUENNE_EPS) ? 1. : fabs (x0/y0 - w0));
  }

  // w and y are bounded (and so is x). Choose between x, y, z
  // depending on intervals first and then to vicinity to midpoint
  CouNumber xl = info -> lower_ [xi], 
            xu = info -> upper_ [xi];

  CouNumber dx = xu-xl,
            dy = yu-yl,
            dw = wu-wl;

  brpts = (double *) realloc (brpts, sizeof (CouNumber));

  // Check largest interval and choose branch variable accordingly.
  // Branching point depends on where the current point is, but for
  // now just focus on the width of the intervals

  way = TWO_RAND;

  if (dx > dy)
    if (dx > dw) {ind = xi; *brpts = (xl + xu) / 2.; return fabs (x0 - y0*w0);} // dx maximum
    else         {ind = wi; *brpts = (wl + wu) / 2.; return fabs (w0 - x0/y0);} // dw 
  else
    if (dy > dw) {ind = yi; *brpts = (yl + yu) / 2.; return fabs (y0 - x0/w0);} // dy
    else         {ind = wi; *brpts = (wl + wu) / 2.; return fabs (w0 - x0/y0);} // dw
}
