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


/// generic approach for negative powers (commom with exprInv::selectBranch
CouNumber negPowSelectBranch (const CouenneObject *obj, 
			      double * &brpts, 
			      int &way,
			      CouNumber k,
			      CouNumber x0, CouNumber y0, 
			      CouNumber l,  CouNumber u) {

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

  // two cases: 
  // 
  // 1) bounds include 0: three way branch to exclude 0 (refer to branchExprDiv.cpp)
  // 2) otherwise         two   way branch

  //  int wi = obj -> Reference () -> Index ();

  if ((l < 0) && (u > 0)) {

    brpts = (double *) realloc (brpts, sizeof (CouNumber));
    *brpts = 0;
    way = TWO_RAND;

    // Closest branch of the hyperbola is on the same side of y+x=0 as
    // (x0,y0) => need only one powNewton

    if (fabs (x0) < COUENNE_EPS)
      x0 = (x0 < 0) ? -COUENNE_EPS : COUENNE_EPS;

    CouNumber xp, xx0 = x0, yy0 = y0, exponent = k;

    if ((x0+y0 < 0) && (x0 > 0) || 
	(x0+y0 > 0) && (x0 < 0)) {

      exponent = 1. / k;
      xx0 = y0;
      yy0 = x0;
    }

    powertriplet pt (exponent);

    xp = (xx0 >= 0) ? 
       powNewton (xx0, yy0, &pt) : 
      -powNewton (-xx0, -yy0, &pt);

    x0 -= xp;
    y0 -= safe_pow (xp, 1. / k);

    return sqrt (x0*x0 + y0*y0); // exact distance
  }

  // case 2: look if inside or outside of belly (refer to branchExprExp.cpp)

  if ((x0 >= 0) && (y0 <  safe_pow  (x0,k)) ||
      (x0 <  0) && (y0 > -safe_pow (-x0,k))) { // outside bellies 

    way = (u < 0) ? TWO_RIGHT : TWO_LEFT; // explore finite interval first

    powertriplet pt (k); // TODO: there may (will) be a problem with negative x0
    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = midInterval ((x0 >= 0) ? 
			   powNewton ( x0,  y0, &pt) : 
			  -powNewton (-x0, -y0, &pt), l, u);

    CouNumber dy = y0 - safe_pow (*brpts >= 0 ? *brpts : - *brpts, 1. / k);
    x0 -= *brpts;
    return sqrt (x0*x0 + dy*dy); // distance is exact
  }

  // Inside, x0*y0 >= 1. Two cases: /////////////////////////////////////////////////
 
  // 1) bounds are infinite both horizontally and vertically
  //    (i.e. [-inf,0] or [0,inf]) --> three way branching

  if ((l <   COUENNE_EPS) && (u >   COUENNE_INFINITY) || 
      (u > - COUENNE_EPS) && (l < - COUENNE_INFINITY)) {

    brpts = (double *) realloc (brpts, 2 * sizeof (double));
    way = THREE_CENTER; // focus on central convexification first

    brpts [0] = x0;      // draw vertical   from (x0,y0) south (north) to curve y=1/x
    brpts [1] = 1. / y0; //      horizontal              west  (east)

    CouNumber a = fabs (y0 - 1 / x0), // sides of a triangle with (x0,y0)
              b = fabs (x0 - 1 / y0), // as one of the vertices
              c = a * cos (atan (a/b));

    return CoinMin (a, CoinMin (b, c)); // distance is exact
  }

  // 2) at least one of them is finite --> two way branching

  brpts = (double *) realloc (brpts, sizeof (double));

  if (l < - COUENNE_INFINITY) { // u << -0

    way = TWO_RIGHT;
    *brpts = midInterval (x0, l, u);

    return CoinMin (y0 - safe_pow (*brpts, 1. / k), 
		    projectSeg (x0, y0, l, safe_pow (l, k), 
				*brpts, safe_pow (*brpts, k), -1)); // distance is exact
  }

  if (u > COUENNE_INFINITY) { // l >> +0

    way = TWO_LEFT;
    *brpts = midInterval (x0, l, u);

    return CoinMin (y0 - safe_pow (*brpts, 1. / k), 
		    projectSeg (x0, y0, l, safe_pow (l, k), 
				*brpts, safe_pow (*brpts, k), +1)); // distance is exact
  }

  // last case: nice finite interval and limited curve

  powertriplet ft (k);
  *brpts = obj -> getBrPoint (&ft, x0, y0, l, u);

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

  return sqrt (x0*x0 + y0*y0); // distance is exact
}



/// set up branching object by evaluating many branching points for
/// each expression's arguments

CouNumber exprInv::selectBranch (const CouenneObject *obj, 
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts, 
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

  return negPowSelectBranch (obj, brpts, way, -1, x0, y0, l,  u);
}
