/*
 * Name:    branchExprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for powers
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprPow.hpp>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprPow::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {

  // return branching point and branching policies for an expression
  // of the form x^k

  ind    = arglist_ [0] -> Index ();
  int wi = w            -> Index ();

  if ((ind < 0) || (wi < 0)) 
    {printf ("Couenne, w=exp(x): negative index\n"); exit (-1);}

  if (arglist_ [1] -> Type () != CONST) {

    printf ("exprPow::selectBranch: warning, power has non-constant exponent\n");
    ind = -1;
    return 0.;
  }

  double k = arglist_ [1] -> Value ();

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  // case 1: k is positive, even, no bounds on x

  int intk = 0;

  bool isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  if (isInt && !(intk % 2)) {

    if ((l < - COUENNE_INFINITY) && 
	(u >   COUENNE_INFINITY)) {

      double alpha = pow ((y0 + pow (x0, k))/2, 1./k),
             yroot = pow (y0, 1./k);

      brpts = (double *) realloc (brpts, 2 * sizeof (double));

      double lambdaL = (-x0 / yroot), lambdaR = 0;

      if (lambdaL < 0) {
	lambdaR = -lambdaL;
	lambdaL = 0;
      }

      brpts [0] = -alpha + lambdaL * (alpha - yroot);
      brpts [1] =  alpha + lambdaR * (yroot - alpha);

      way = THREE_CENTER;
    }
  }

  ind = -1;
  return 0.;

  // the main cases depend on the bounds.

  // case

  //
  // 1) x^k is convex
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
  /*
  ind    = argument_ -> Index ();
  int wi = w         -> Index ();

  if ((ind < 0) || (wi < 0)) {printf ("Couenne, w=exp(x): negative index\n"); exit (-1);}

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  if (y0 < exp (x0)) { 

    // Outside

    brpts = (double *) realloc (brpts, sizeof (double));

    *brpts = powNewton (x0, y0, exp, exp, exp);

    if      (*brpts < l) *brpts = (u <   COUENNE_INFINITY) ? ((l+u)/2) : l + 1;
    else if (*brpts > u) *brpts = (l > - COUENNE_INFINITY) ? ((l+u)/2) : u - 1;

    way = TWO_RAND;
    CouNumber dy = y0 - exp (*brpts);
    x0 -= *brpts;
    return sqrt (x0*x0 + dy*dy);

  } else { 

    // Inside. Two cases:
 
    if ((l < -COUENNE_INFINITY) && 
	(u >  COUENNE_INFINITY)) {

      // 1) both bounds are infinite --> three way branching

      brpts = (double *) realloc (brpts, 2 * sizeof (double));
      way = THREE_CENTER; // focus on central convexification first

      brpts [0] = x0;       // draw vertical   from (x0,y0) south to curve y=exp(x)
      brpts [1] = log (y0); //      horizontal              east

      CouNumber a = y0 - exp (x0), // sides of a triangle with (x0,y0)
		b = x0 - log (y0), // as one of the vertices
        	c = a * cos (atan (a/b));

      return mymin (a, mymin (b, c));

    } else {

      // 2) at least one of them is finite --> two way branching

      brpts = (double *) realloc (brpts, sizeof (double));

      if (l < -COUENNE_INFINITY) {

	*brpts = x0;
	if (*brpts > u - COUENNE_NEAR_BOUND) *brpts = u-1;
	way = TWO_RIGHT;
	return y0 - exp (x0);

      } else if (u > COUENNE_INFINITY) {

	*brpts = log (y0);
	if (*brpts < l + COUENNE_NEAR_BOUND) *brpts = l+1;
	way = TWO_LEFT;
	return x0 - log (y0);

      } else {

	*brpts = powNewton (x0, y0, exp, exp, exp);

	if ((*brpts > u - COUENNE_NEAR_BOUND) ||
	    (*brpts < l + COUENNE_NEAR_BOUND))
	  *brpts = (l+u) / 2;

	way = TWO_RAND;

	CouNumber dx = x0 - *brpts,
  	          dy = y0 - exp (*brpts);

	return sqrt (dx*dx + dy*dy);
      }
    }
  }
  */

  // TODO...

  //  ind = argument_ -> Index ();

  /*  if (ind < 0) {
    printf ("argument of exprAbs has negative index\n");
    exit (-1);
  }

  brpts = (double *) malloc (sizeof (double));
  *brpts = 0.;
  way = TWO_LEFT;

  return mymin (project (1., -1., 0., x0, y0, 0., COUENNE_INFINITY,  0, NULL, NULL),
  project (1., +1., 0., x0, y0, -COUENNE_INFINITY, 0., 0, NULL, NULL));*/
}
