/*
 * Name:    branchExprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for exponentials
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprExp.h>
#include <exprPow.h>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprExp::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {
  ind = -1;
  return 0.;

  // two cases: inside or outside the belly. 
  //
  // Inside: the distance depends on the projection of the current
  // point onto the would-be upper envelopes, which forces us to look
  // at it numerically. If both bounds are infinite, create a ThreeWay
  // branch.
  //
  // Outside: it suffices to project the current point on the line to
  // get the maxi-min displacement.
  //
  // As happens for all monotonous functions, after chosing *brpts it
  // is equivalent to choose w's or x's index as ind, as the implied
  // bounds will do the rest.

  ind = argument_ -> Index ();
  int wi = w -> Index ();

  if ((ind < 0) || (wi < 0)) {printf ("Couenne, w=f(x): negative index\n"); exit (-1);}

  brpts = (double *) realloc (brpts, sizeof (double));

  CouNumber x0 = info -> solution_ [ind],
            y0 = info -> solution_ [wi];

  if (y0 < exp (x0)) { // A) Outside

    *brpts = powNewton (x0, y0, exp, exp, exp);
    way = TWO_RAND;
    CouNumber dy = y0 - exp (*brpts);
    x0 -= *brpts;
    return sqrt (x0*x0 + dy*dy);

  } else { // B) Inside 

    CouNumber l = info -> lower_ [ind],
              u = info -> upper_ [ind];

    // two cases:
 
    if ((l < -COUENNE_INFINITY) && 
	(u >  COUENNE_INFINITY)) {

      // B1) bounds are infinite

      brpts = (double *) realloc (brpts, 2*sizeof (double));
      way = THREE_CENTER;
      *brpts = x0;
      brpts [1] = log (y0);

      CouNumber a = y0 - exp (x0), // sides of a triangle with (x0,y0)
		b = x0 - log (y0), // as one of the vertices
        	c = a * cos (atan (a/b));

      return mymin (a, mymin (b, c));

    } else {

      brpts = (double *) realloc (brpts, sizeof (double));

      // B2) at least one of them is finite

      if (l < -COUENNE_INFINITY) {

	*brpts = x0;
	way = TWO_RIGHT;
	return y0 - exp (x0);

      } else if (u > COUENNE_INFINITY) {

	*brpts = log (y0);
	way = TWO_LEFT;
	return x0 - log (y0);

      } else {

	*brpts = powNewton (x0, y0, exp, exp, exp);
	way = TWO_RAND;

	CouNumber dx = x0 - *brpts,
  	          dy = y0 - exp (*brpts);

	return sqrt (dx*dx + dy*dy);
      }
    }
  }
}
