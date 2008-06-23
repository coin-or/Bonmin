/*
 * Name:    branchExprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch selection for 1/x
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprInv.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "projections.hpp"
#include "funtriplets.hpp"


/// generic approach for negative powers (commom with exprInv::selectBranch())
CouNumber negPowSelectBranch (const CouenneObject *obj, 
			      double * &brpts, 	
			      double * &brDist, // distance of current LP
			                        // point to new convexifications
			      int &way,
			      CouNumber k,
			      CouNumber x0, CouNumber y0, 
			      CouNumber l,  CouNumber u) {

  brDist = (double *) realloc (brDist, 2 * sizeof (double));
  brpts  = (double *) realloc (brpts, sizeof (CouNumber));

  // two cases: inside or outside the curves (there are two branches
  // of the hyperbola).
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

  if ((l < -COUENNE_EPS) && (u > COUENNE_EPS)) { // handle discontinuity

    // no matter if the (negative) exponent is odd or even, we better
    // branch on 0 lest we have no good convexification (especially
    // with odd exponent)

    *brpts = 0.;
    way = TWO_RAND;

    // Closest branch of the hyperbola is on the same side of y+x=0 as
    // (x0,y0) => need only one powNewton

    if (fabs (x0) < COUENNE_EPS)
      x0 = (x0 <= -0.) ? -COUENNE_EPS : COUENNE_EPS;

    CouNumber xp, xx0 = x0, yy0 = y0, exponent = k;

    // invert dependent and independent if
    if ((x0+y0 < 0.) && (x0 > 0.) ||  // in lower half of fourth orthant, or
	(x0+y0 > 0.) && (x0 < 0.)) {  // in upper half of second orthant

      exponent = 1. / k;
      xx0 = y0;
      yy0 = x0;
    }

    powertriplet pt (exponent);

    xp = (xx0 >= 0) ? 
       powNewton  (xx0,  yy0, &pt) : 
      -powNewton (-xx0, -yy0, &pt);

    CouNumber diff = x0 - xp;
    y0 -= safe_pow (xp, 1. / k);

    brDist [0] = sqrt (diff*diff + y0*y0); // exact distance
    brDist [1] = CoinMax (fabs (x0), 1.);

    if (x0 > 0.) {
      double swap = brDist [0];
      brDist [0] = brDist [1];
      brDist [1] = swap;
    }

    return CoinMin (brDist [0], brDist [1]);
  }

  int intk = 0;

  bool
    isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
    isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  // case 2: bound interval does not contain zero. Look if inside or
  // outside of belly (refer to branchExprExp.cpp)

  if ((x0 >=  0.) &&                          (y0 <  safe_pow  (x0,k))  ||   // x0>0, or
      (x0 <= -0.) &&                                                         // x0<0, and
      ((isInt &&               !(intk % 2) && (y0 <  safe_pow  (x0,k))) ||     // integer, even
       ((isInt || isInvInt) &&  (intk % 2) && (y0 > -safe_pow (-x0,k))))) {    // (inv)integer, odd

    // Inside. Branch on closest point on curve, computed with a
    // Newton method

    way = (u < -0.) ? TWO_RIGHT : TWO_LEFT; // explore finite interval first

    powertriplet pt (k);

    *brpts = obj -> midInterval ((x0 >= 0.) ? 
 	 			  powNewton ( x0,  y0, &pt) : 
				 -powNewton (-x0, -y0, &pt), l, u);

    CouNumber dy = y0 - safe_pow (*brpts >= 0 ? *brpts : - *brpts, 1. / k);
    x0 -= *brpts;
    return (brDist [0] = brDist [1] = sqrt (x0*x0 + dy*dy)); // distance is exact
  }

  // Inside, (x0^k) * y0 >= 1. Two cases: /////////////////////////////////////////////////
 
  // 1) bounds are infinite both horizontally and vertically
  // (i.e. [-inf,0] or [0,inf]) --> as for exprExp, pick point on
  // diagonal from current to curve, to be sure current will be cut by
  // branching rule

  if ((l <   COUENNE_EPS) && (u >   COUENNE_INFINITY) || 
      (u > - COUENNE_EPS) && (l < - COUENNE_INFINITY)) {

    /* brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first
    brpts [0] = x0;      // draw vertical   from (x0,y0) south (north) to curve y=1/x
    brpts [1] = 1. / y0; //      horizontal              west  (east)
    CouNumber a = fabs (y0 - 1 / x0), // sides of a triangle with (x0,y0)
              b = fabs (x0 - 1 / y0), // as one of the vertices
              c = a * cos (atan (a/b)); */

    //brpts = (double *) realloc (brpts, sizeof (double));

    //if (x0 > COUENNE_EPS) 
    *brpts = 0.5 * (fabs (x0) + pow (fabs (y0), 1./k));

    if (x0 < 0.) {
      *brpts = - *brpts;
      brDist [0] = fabs (fabs (y0) - safe_pow (fabs (x0), k));
      brDist [1] = *brpts - x0;
    } else {
      brDist [0] = x0 - *brpts;
      brDist [1] = fabs (y0 - safe_pow (x0, k));
    }

    //else 
    //*brpts = 0.5 * (x0 + pow (y0, 1./k));

    // follow South-East diagonal to find point on curve
    // so that current point is surely cut 
    //*brpts = 0.5 * (x0 + log (y0)); 
    //way = TWO_RAND;
    way = (x0 > *brpts) ? TWO_RIGHT : TWO_LEFT;

    return CoinMin (brDist [0], brDist [1]);
    //x0 - pow (fabs (y0), 1./k), y0 - pow (x0,k));
    //return CoinMin (a, CoinMin (b, c)); // distance is exact
  }

  // 2) at least one of them is finite

  if (l < - COUENNE_INFINITY) { // u << -0

    way = TWO_RIGHT;
    *brpts = obj -> midInterval (x0, l, u);

    return CoinMin (brDist [0] = y0 - safe_pow (*brpts, 1. / k), 
		    brDist [1] = projectSeg (x0, y0, l, safe_pow (l, k), 
					     *brpts, safe_pow (*brpts, k), -1)); // distance is exact
  }

  if (u > COUENNE_INFINITY) { // l >> +0

    way = TWO_LEFT;
    *brpts = obj -> midInterval (x0, l, u);

    return CoinMin (brDist [1] = y0 - safe_pow (*brpts, 1. / k), 
		    brDist [0] = projectSeg (x0, y0, l, safe_pow (l, k), 
					     *brpts, safe_pow (*brpts, k), +1)); // distance is exact
  }

  // last case: nice finite interval and limited curve

  powertriplet ft (k);
  *brpts = obj -> getBrPoint (&ft, x0, l, u);

  /*  // TODO: check if it works with all exponents
  if (u > l + COUENNE_EPS) {

    powertriplet ft (k);
    *brpts = maxHeight (&ft, l, u); // min area

    // *brpts = safe_pow ((safe_pow (u,k) - safe_pow (l,k)) / (k * (u-l)), 1/(k-1));
    // if (u < 0)
    // *brpts = - *brpts;
  }
  else *brpts = midInterval (x0, l, u);*/

  way = TWO_RAND;

  x0 -= *brpts;
  y0 -= safe_pow (*brpts, k);

  brDist [0] = projectSeg (x0,y0, l,      safe_pow (l,      k), *brpts, safe_pow (*brpts, k), 0);
  brDist [1] = projectSeg (x0,y0, *brpts, safe_pow (*brpts, k), u,      safe_pow (u,      k), 0);

  return CoinMin (brDist [0], brDist [1]);//sqrt (x0*x0 + y0*y0); // distance is exact
}



/// set up branching object by evaluating many branching points for
/// each expression's arguments

CouNumber exprInv::selectBranch (const CouenneObject *obj, 
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts, 
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {

  var = argument_;

  int
    ind = argument_           -> Index (),
    wi  = obj -> Reference () -> Index ();

  assert ((ind >= 0) && (wi >= 0));

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  return negPowSelectBranch (obj, brpts, brDist, way, -1, x0, y0, l,  u);
}
