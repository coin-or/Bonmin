/*
 * Name:    branchExprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for multiplications
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprMul.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "domain.hpp"

#define LARGE_BOUND 9.9e12
#define BOUND_WING 1

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprMul::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts, 
				 int &way) {

  Domain *domain = arglist_ [0] -> domain ();

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index (),
      wi = obj -> Reference () -> Index ();

  assert ((xi >= 0) && (yi >= 0) && (wi >= 0));

  CouNumber 
    x0 = domain -> x  (xi), y0 = domain -> x  (yi),
    xl = domain -> lb (xi), yl = domain -> lb (yi),
    xu = domain -> ub (xi), yu = domain -> ub (yi),
    w0 = domain -> x (wi);

  // First, try to avoid infinite bounds for multiplications, which
  // make them pretty hard to deal with

  if (((var = arglist_ [0]) -> Index() >= 0) && (xl < -COUENNE_INFINITY) && (xu > COUENNE_INFINITY) ||
      ((var = arglist_ [1]) -> Index() >= 0) && (yl < -COUENNE_INFINITY) && (yu > COUENNE_INFINITY)) {

    // branch around current point. If it is also at a crazy value,
    // reset it close to zero.

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    CouNumber curr = (*var) ();//expression::Variable (ind);

    if (fabs (curr) >= LARGE_BOUND) curr = 0;

    //brpts [0] = - fabs (curr) - BOUND_WING;
    //brpts [1] =   fabs (curr) + BOUND_WING;
    brpts [0] = curr - BOUND_WING;
    brpts [1] = curr + BOUND_WING;

    way = THREE_CENTER;

    return fabs (w0 - x0*y0);
    //    return - COUENNE_INFINITY; // tell caller not to set infeasibility to this
  }

  // at least one bound is infinite

  brpts = (double *) realloc (brpts, sizeof (double)); // only one branch point

  // don't privilege xi over yi

  int ind;

  if      ((xl < -COUENNE_INFINITY) && (xu > COUENNE_INFINITY)) // x unbounded in both directions
    {ind = xi; *brpts = 0; way = TWO_RAND;}

  else if ((yl < -COUENNE_INFINITY) && (yu > COUENNE_INFINITY)) // y unbounded in both directions
    {ind = yi; *brpts = 0; way = TWO_RAND;}

  else if (xl < -COUENNE_INFINITY)                              // x unbounded below
    {ind = xi; *brpts = obj -> midInterval (((x0 < 0.) ? 2 : 0.5) * x0, xl, xu); way = TWO_RIGHT;}

  else if (xu >  COUENNE_INFINITY)                              // x unbounded above
    {ind = xi; *brpts = obj -> midInterval (((x0 > 0.) ? 2 : 0.5) * x0, xl, xu); way = TWO_LEFT;} 

  else if (yl < -COUENNE_INFINITY)                              // y unbounded below
    {ind = yi; *brpts = obj -> midInterval (((y0 < 0.) ? 2 : 0.5) * y0, yl, yu); way = TWO_RIGHT;}

  else if (yu >  COUENNE_INFINITY)                              // y unbounded above
    {ind = yi; *brpts = obj -> midInterval (((y0 > 0.) ? 2 : 0.5) * y0, yl, yu) ;way = TWO_LEFT;} 

  else { // both are bounded

    //way = TWO_RAND;

    CouNumber delta = (yu-yl) - (xu-xl);

    if      (delta > +COUENNE_EPS) ind = yi;
    else if (delta < -COUENNE_EPS) ind = xi;
    else ind = (CoinDrand48 () < 0.5) ? xi : yi;

    CouNumber 
      pt = domain -> x  (ind),
      lb = domain -> lb (ind),
      ub = domain -> ub (ind);

#define THRES_ZERO_SYMM 0.7

    if ((lb < -COUENNE_EPS) && (ub > COUENNE_EPS) && 
	(-lb/ub >= THRES_ZERO_SYMM) &&
	(-ub/lb >= THRES_ZERO_SYMM))
      // interval is fairly symmetric around 0, branch on it
      *brpts = 0.;
    else
      *brpts = (obj -> Strategy () != CouenneObject::MID_INTERVAL) ? 
	(0.5 * (lb+ub)) :
	obj -> midInterval (pt, lb, ub);

    way = (pt > *brpts) ? TWO_RIGHT : TWO_LEFT;
  }

  var = arglist_ [(ind == xi) ? 0 : 1];

  return fabs (w0 - x0*y0);
}
