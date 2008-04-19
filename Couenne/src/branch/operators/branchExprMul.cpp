/*
 * Name:    branchExprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for multiplications
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

#include "exprMul.hpp"
#include "funtriplets.hpp"
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

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index (),
      wi = obj -> Reference () -> Index ();

  assert ((xi >= 0) && (yi >= 0) && (wi >= 0));

  CouNumber 
    x0 = info -> solution_  [xi], y0 = info -> solution_  [yi],
    xl = info -> lower_     [xi], yl = info -> lower_     [yi],
    xu = info -> upper_     [xi], yu = info -> upper_     [yi],
    w0 = info -> solution_  [wi];

#ifdef DEBUG
  printf ("    branch MUL: %g [%g,%g] %g [%g,%g]\n", 
	  x0, xl, xu, y0, yl, yu);
#endif

  brpts = (double *) realloc (brpts, sizeof (double));

  if (fabs (xu-xl) < COUENNE_EPS) { // x almost constant

    if (fabs (yu-yl) < COUENNE_EPS) { // both almost constant, return null result

      var = NULL;
      return 0.;

    } else { // x only, branch on y

      var = arglist_ [1];
      *brpts = 0.5 * (yl+yu);
      return fabs (w0 - x0*y0);
    }

  } else if (fabs (yu-yl) < COUENNE_EPS) { // y almost constant, branch on x

    var = arglist_ [0];
    *brpts = 0.5 * (xl+xu);
    return fabs (w0 - x0*y0);
  }

  // First, try to avoid infinite bounds for multiplications, which
  // make them pretty hard to deal with

  if (((var = arglist_ [0]) -> Index() >= 0) && (xl < -COUENNE_INFINITY) && (xu > COUENNE_INFINITY) ||
      ((var = arglist_ [1]) -> Index() >= 0) && (yl < -COUENNE_INFINITY) && (yu > COUENNE_INFINITY)) {

    // restore it when we have three-way branching

#if 0
    // branch around current point. If it is also at a crazy value,
    // reset it close to zero.

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    CouNumber curr = (*var) ();

    if (fabs (curr) >= LARGE_BOUND) curr = 0;

    //brpts [0] = - fabs (curr) - BOUND_WING;
    //brpts [1] =   fabs (curr) + BOUND_WING;
    brpts [0] = curr - BOUND_WING;
    brpts [1] = curr + BOUND_WING;

    way = THREE_CENTER;
#endif

    *brpts = 0.;

    return fabs (w0 - x0*y0);
  }

  // don't privilege xi over yi

  int ind;

  // at least one bound is infinite

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
      pt = info -> solution_  [ind],
      lb = info -> lower_     [ind],
      ub = info -> upper_     [ind];

    if ((lb < -COUENNE_EPS) && (ub > COUENNE_EPS) && 
	(-lb/ub >= THRES_ZERO_SYMM) &&
	(-ub/lb >= THRES_ZERO_SYMM))
      // interval is fairly symmetric around 0, branch on it
      *brpts = 0.;

    else switch (obj -> Strategy ()) {
      case CouenneObject::MID_INTERVAL: *brpts = obj -> midInterval (pt, lb, ub);             break;
      case CouenneObject::BALANCED:     *brpts = balancedMul (info, (ind == xi) ? 0 : 1, wi); break;
      case CouenneObject::MIN_AREA:
      default:                          *brpts = (0.5 * (lb+ub));                             break;
    }

    way = (pt > *brpts) ? TWO_RIGHT : TWO_LEFT;
  }

  var = arglist_ [(ind == xi) ? 0 : 1];

#ifdef DEBUG
  printf ("    MUL: br on x_%d %g [%g,%g] [%g,%g] (%g,%g)\n", 
	  ind, *brpts, xl, xu, yl, yu, x0, y0);
#endif

  return fabs (w0 - x0*y0);
}


// branching point for multiplication according to the balanced strategy
CouNumber exprMul::balancedMul (const OsiBranchingInformation *info, int index, int wind) {

  // first of all, make sure both arguments are variables

  int other;

  if (index==0) {
    index = arglist_ [0] -> Index ();
    other = arglist_ [1] -> Index ();
  } else {
    index = arglist_ [1] -> Index ();
    other = arglist_ [0] -> Index ();
  }

  assert ((index >= 0) && (other >= 0));

  CouNumber 
    xl = info -> lower_    [index],  yl = info -> lower_    [other],
    xu = info -> upper_    [index],  yu = info -> upper_    [other],
    x0 = info -> solution_ [index],  y0 = info -> solution_ [other],
    w0 = info -> solution_ [wind];
    
  // It is quite tricky to implement a balanced strategy for products,
  // because it is more difficult to measure "balancedness" for binary
  // operators than it is for univariate functions.
  //
  // As a rule of thumb, we therefore apply the usual balanced
  // strategy for the univariate function resulting from constraining
  // (x,y) to a segment crossing their bounding box [xl,xu] X [yl,yu].
  //
  // Said segment is the set of points between (xl,yl) and (xu,yu) if
  // the current point is above the curve w:=xy, otherwise it is the
  // other diagonal, i.e. the set of points between (xl,yu) and
  // (xu,yl).
  //
  // In the two cases, we have the point 
  //
  // (above) P(t) = (xp,yp) := (xl + t (xu-xl), yl + t (yu-yl))
  // (below) P(t) = (xp,yp) := (xu + t (xl-xu), yl + t (yu-yl))
  //
  // with t in [0,1], which forms the following second degree
  // polynomial when multiplying the coordinates:
  //
  // (above) f(t) = xp*yp = (yu-yl)*(xu-xl)*t^2 + [yl(xu-xl) + xl(yu-yl)]*t + xl*yl
  // (below) f(t) = xp*yp = (yu-yl)*(xl-xu)*t^2 + [yl(xl-xu) + xu(yu-yl)]*t + xu*yl
  //
  // which is a quadratic function that can be expressed in the form
  // of f'(z) = z^2 + c if we apply an affine transformation to t:
  //
  // t = mz + q
  //
  // such that the resulting coefficients of the quadratic- and the
  // linear terms are one and zero, respectively. Thus:
  // 
  // (above) f(z) = (yu-yl)*(xu-xl)*(mz+q)^2               + [yl(xu-xl)-xl(yu-yl)]*(mz+q) + xl*yl = 
  //              = (yu-yl)*(xu-xl)*(m^2 z^2 + 2qmz + q^2) + [yl(xu-xl)-xl(yu-yl)]*(mz+q) + xl*yl = 
  //              = z^2 + c 
  //
  // (below) f(z) = (yu-yl)*(xl-xu)*(mz+q)^2               + [yl(xl-xu)-xu(yu-yl)]*(mz+q) + xu*yl = 
  //              = (yu-yl)*(xl-xu)*(m^2 z^2 + 2qmz + q^2) + [yl(xl-xu)-xu(yu-yl)]*(mz+q) + xu*yl = 
  //              = -z^2 + c 
  //
  // if and only if
  //
  // (above) ((yu-yl)*(xu-xl)) * m^2 =  1   )
  //                                        } <====>  m = 1 / sqrt ((yu-yl)*(xu-xl))
  // (below) ((yu-yl)*(xl-xu)) * m^2 = -1   )
  //
  // (above)  2qm*(yu-yl)*(xu-xl) + m[yl(xu-xl)-xl(yu-yl)] = 0   <====>
  //          q = -[yl(xu-xl)-xl(yu-yl)] / (2(yu-yl)*(xu-xl))
  //
  // (below)  2qm*(yu-yl)*(xl-xu) + m[yl(xl-xu)-xu(yu-yl)] = 0   <====>
  //          q = -[yl(xl-xu)-xu(yu-yl)] / (2(yu-yl)*(xl-xu))
  //
  // If the point is below the curve, a very similar reasoning applies
  // (simply swap xl with xu).
  // 
  // Hence, we simply apply the balanced strategy to the function
  // f(z)=z^2 with z bounded between -q/m and (1-q)/m. The returning
  // value z_opt must be transformed to get
  //
  // t_opt = m z_opt + q
  //
  // and the branching point is xl + t_opt (xu-xl)
  //

  powertriplet ft (2);

  // A: above
  // B: below

  bool above = (w0 > x0*y0);

  CouNumber 
    dx     = xu-xl,
    dy     = yu-yl,
    area   = dx*dy,          // coefficient of t^2
    bA     =  yl*dx - xl*dy, // coefficient of t
    bB     = -yl*dx - xu*dy, // coefficient of t
    m      = 1. / sqrt (area),
    qA     = -bA / (2*area),
    qB     =  bB / (2*area),
    z_opt = above ? 
      minMaxDelta (&ft, -qA/m, (1-qA)/m): 
      minMaxDelta (&ft, -qB/m, (1-qB)/m),
    t_optA = m*z_opt + qA,
    t_optB = m*z_opt + qB;

  /*printf ("------------------\n(%d,%d): [%g,%g], [%g,%g]\n", index, other, xl, xu, yl, yu);
  printf ("dx = %g, dy = %g, area = %g, bA = %g, bB = %g, m = %g\n", dx, dy, area, bA, bB, m);
  printf ("qA = %g, qB = %g, z_opt = %g, tA = %g, tB = %g\n", 
	  qA, qB, z_opt, t_optA, t_optB);
  printf ("w = %g %c xy = %g*%g = %g --> %g\n", w0,
	  (w0 > x0*y0) ? '>' : '<', x0, y0, x0*y0, 
	  (w0 > x0*y0) ? (xl + dx * t_optA) : (xu - dx * t_optB));*/

  return (w0 > x0*y0) ? 
    (xl + dx * t_optA) : 
    (xu - dx * t_optB);
}
