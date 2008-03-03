/*
 * Name:    branchExprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for logarithms
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprLog.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "projections.hpp"
#include "funtriplets.hpp"


#define SQ_COUENNE_EPS COUENNE_EPS * COUENNE_EPS

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprLog::selectBranch (const CouenneObject *obj, 
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts, 
				 int &way) {

  // quite similar to exprExp::selectBranch() (see branchExprExp.cpp)
  //
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

  var = argument_;

  int
    ind = var -> Index (),
    wi  = obj -> Reference () -> Index ();

  assert ((ind >= 0) && (wi >= 0));

  CouNumber 
    y0 = info -> solution_ [wi],
    x0 = info -> solution_ [ind],
    l  = info -> lower_    [ind],
    u  = info -> upper_    [ind];

  if (u < COUENNE_EPS) { // strange case, return default branching rule
    var = NULL; 
    return 0.;
  }

  if (x0 < SQ_COUENNE_EPS) // very unlikely...
    x0 = SQ_COUENNE_EPS;

  if (y0 > log (x0)) { 

    // Outside -> branch on closest point on curve

    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = obj -> midInterval (powNewton (x0, y0, log, inv, oppInvSqr), l, u);

    way = TWO_LEFT;
    CouNumber dy = y0 - log (*brpts);
    x0 -= *brpts;

    return sqrt (x0*x0 + dy*dy); // exact distance
  } 

  // Inside. Two cases:
 
  if ((l <= SQ_COUENNE_EPS) && 
      (u > COUENNE_INFINITY)) {

    // 1) curve is unlimited in both senses --> three way branching

    // restore when we can do three-way branching

#if 0
    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = exp (y0); // draw horizontal from (x0,y0) east  to curve y=log(x)
    brpts [1] = x0;       //      vertical                north

    CouNumber a = x0 - exp (y0), // sides of a triangle with (x0,y0)
              b = log (x0) - y0; // as one of the vertices

    return a * cos (atan (a/b)); // exact distance
#endif

    // follow South-East diagonal to find point on curve
    // so that current point is surely cut 
    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = 0.5 * (x0 + exp (y0));
    way = TWO_RAND;

    return CoinMin (fabs (x0 - exp(y0)), fabs (y0 - log (x0)));
  } 

  // 2) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l <= SQ_COUENNE_EPS) { // u is finite

    *brpts = obj -> midInterval (exp (y0), l, u);
    way = TWO_RIGHT;

  return projectSeg (x0, y0, *brpts, log (*brpts), u, log (u), +1); // exact distance
    //    return CoinMin (x0 - exp (y0), log (x0) - y0);
  }
 
  if (u > COUENNE_INFINITY) { // l is far from zero

    *brpts = obj -> midInterval (x0, l, u);
    way = TWO_LEFT;

    return projectSeg (x0, y0, l, log (l), *brpts, log (*brpts), +1); // exact distance
    //return log (x0) - y0;
  } 

  // both are finite
 
  simpletriplet ft (log, inv, oppInvSqr, inv);

  *brpts = obj -> getBrPoint (&ft, x0, l, u);

  //  *brpts = midInterval (powNewton (x0, y0, log, inv, oppInvSqr), l, u); 
  // WRONG! Local minima may be at bounds

  // compute distance from current point to new convexification(s) and
  // to curve. If closer to curve, branch on current point

  way = TWO_RAND;

  // exact distance
  return CoinMin (projectSeg (x0, y0, l, log (l), *brpts, log (*brpts),             +1),
		  projectSeg (x0, y0,             *brpts, log (*brpts), u, log (u), +1));
}
