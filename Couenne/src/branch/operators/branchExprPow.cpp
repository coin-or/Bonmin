/*
 * Name:    branchExprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for powers
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprPow.hpp"
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
			      CouNumber l,  CouNumber u);


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprPow::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression *&var,
				 double * &brpts, 
				 int &way) {

  // return branching point and branching policies for an expression
  // of the form x^k

  var = arglist_ [0];

  int 
    ind = arglist_ [0]        -> Index (),
    wi  = obj -> Reference () -> Index ();

  assert ((ind >= 0) && (wi >= 0) && (arglist_ [1] -> Type () == CONST));

  double k = arglist_ [1] -> Value ();

  CouNumber y0 = info -> solution_ [wi],
            x0 = info -> solution_ [ind],
            l  = info -> lower_    [ind],
            u  = info -> upper_    [ind];

  /*printf ("selbra for "); print ();
  printf ("%g [%g,%g] -> %g [%g,%g]\n", 
	  x0, l, u, y0, 
	  info -> lower_    [wi],
	  info -> upper_    [wi]);*/

  if      (x0 < l) x0 = l;
  else if (x0 > u) x0 = u;

  // bounds coincinde (happens within setupList)
  if (fabs (u-l) < COUENNE_EPS) {
    brpts = (CouNumber *) realloc (brpts, sizeof (CouNumber));
    way = TWO_RAND;
    return (y0 - pow (x0, k));
  }

  // case 1: k negative, resort to method similar to exprInv:: ///////////////////////////////
  if (k<0) return negPowSelectBranch (obj, brpts, way, k, x0, y0, l, u);

  int intk = 0;

  bool isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  // case 2: k is positive and even /////////////////////////////////////////////////////////

  if (isInt && !(intk % 2)) {

    if ((l < - COUENNE_INFINITY) && 
	(u >   COUENNE_INFINITY)) { // infinite bounds

      if (y0 < pow (x0, k)) { // good side

	brpts = (double *) realloc (brpts, sizeof (double));
	powertriplet pt (k);
	*brpts = powNewton (x0, y0, &pt);

	way = (x0 > 0) ? TWO_RIGHT : TWO_LEFT; // priority to same side as x0

	x0 -= *brpts;
	y0 -= pow (*brpts, k);

	return sqrt (x0*x0 + y0*y0); // exact distance
      }

      // on the bad side

      double xp = pow (y0, 1./k);

      brpts = (double *) realloc (brpts, 2 * sizeof (double));
      brpts [0] = 0.5 * (x0 - xp);
      brpts [1] = 0.5 * (x0 + xp);

      way = THREE_CENTER;

      return CoinMin (x0 - brpts [0], 
		      CoinMin (brpts [1] - x0, 
			       projectSeg (x0, y0, 
					   brpts [0], pow (brpts [0], k),
					   brpts [1], pow (brpts [1], k), -1)));
      // no bounds on x
      /*      double alpha = pow ((y0 + pow (x0, k))/2, 1./k),
             yroot = pow (y0, 1./k);
      brpts = (double *) realloc (brpts, 2 * sizeof (double));
      double lambdaL = (-x0 / yroot), lambdaR = 0;
      if (lambdaL < 0) {
	lambdaR = -lambdaL;
	lambdaL = 0;
      }
      CouNumber // approx distance
	b0 = brpts [0] = -alpha + lambdaL * (alpha - yroot),
	b1 = brpts [1] =  alpha + lambdaR * (yroot - alpha);
      way = THREE_CENTER;
      return CoinMin (projectSeg (x0, y0, b0, pow (b0, k), b1, pow (b1, k), -1), 
		      CoinMin (x0 - b0, b1 - x0));      */
    }

    // at least one bound is finite //////////////////////////////////////////////

    brpts = (double *) realloc (brpts, sizeof (double));

    if (l < -COUENNE_INFINITY) {

      *brpts = - safe_pow (y0, 1. / k);
      way = TWO_RIGHT;
      return CoinMin (x0 - *brpts, 
		      projectSeg (x0,     y0, 
				  *brpts, safe_pow (*brpts, k), 
				  u,      safe_pow (u, k), -1));
    }

    if (u >  COUENNE_INFINITY) {

      *brpts = safe_pow (y0, 1. / k);
      way = TWO_LEFT;
      return CoinMin (*brpts - x0,
		      projectSeg (x0,     y0, 
				  l,      safe_pow (l, k),
				  *brpts, safe_pow (*brpts, k), -1));
    } 

    powertriplet ft (k);
    //*brpts = maxHeight (&ft, x0, y0, l, u);
    *brpts = obj -> getBrPoint (&ft, x0, y0, l, u);

    way = (x0 < *brpts) ? TWO_LEFT : TWO_RIGHT;

    //    w -> print (); printf (" := "); dynamic_cast <exprAux *> (w) -> Image () -> print ();
    /*print ();
    printf (" (%g,%g) [%g,%g] --> %g, distances = %g,%g\n",
	    x0, y0, l, u, *brpts,
	    projectSeg (x0, y0, l,      safe_pow (l,k), *brpts, safe_pow (*brpts,k), 0),
	    projectSeg (x0, y0, *brpts, safe_pow (*brpts,k),      u, safe_pow (u,k), 0));*/

    // Min area  -- exact distance
    CouNumber retval = 
      CoinMin (projectSeg (x0, y0, l,      safe_pow (l,k), *brpts, safe_pow (*brpts,k), 0),
	       projectSeg (x0, y0, *brpts, safe_pow (*brpts,k),      u, safe_pow (u,k), 0));

    //printf ("returning %.10e\n", retval);

    return retval;

    /*      brpts = (double *) realloc (brpts, sizeof (double));
     *brpts = midInterval (x0, l, u);
     way = 
     (l < - COUENNE_INFINITY) ? TWO_RIGHT : 
     (u >   COUENNE_INFINITY) ? TWO_LEFT  : TWO_RAND;

     return fabs (y0 - pow (x0, k)); // not an exact measure
    */
  }

  // from here on, we use two-way branch

  brpts = (double *) realloc (brpts, sizeof (double));
  CouNumber pow0 = pow (x0, k);

  // case 3: k>1 and odd ////////////////////////////////////////////////////////////

  if (isInt && (intk % 2)) {

    way = (x0 > 0) ? TWO_RIGHT : TWO_LEFT;

    if ((l < - COUENNE_INFINITY) && (u > COUENNE_INFINITY) || // [-inf,+inf[
	(l < - COUENNE_INFINITY) && (y0 < pow0)            ||
	(u >   COUENNE_INFINITY) && (y0 > pow0)){ 

	if ((y0 > 0) && (y0 < pow0) ||  
	    (y0 < 0) && (y0 > pow0)) {

	  *brpts = 0;
	  return fabs (pow0 - y0);

	} else {

	  *brpts = pow (y0, 1./k);

	  return (y0 > 0) ? // approx distance
	    projectSeg (x0, y0, x0, CoinMax (pow0, 0.), *brpts, y0, 0) :
	    projectSeg (x0, y0, x0, CoinMin (pow0, 0.), *brpts, y0, 0);
	}
    }

    // otherwise, the convexification is surely bounded. 
    //
    // apply minarea

    if (u > l + COUENNE_EPS) {
      *brpts = safe_pow ((safe_pow (u, k) - safe_pow (l,k)) / (k* (u-l)), 1./(k-1.));
      if (u<0) 
	*brpts = -*brpts;
    }
  }


  // case 4: k positive, in ]0,1[ and 1/k is integer and odd ////////////////////////

  if (isInvInt && (intk % 2)) {

    way = (x0 > 0) ? TWO_RIGHT : TWO_LEFT;

    if ((l < - COUENNE_INFINITY) && (u > COUENNE_INFINITY) || // ]-inf,+inf[
	(l < - COUENNE_INFINITY) && (y0 < pow0)            ||
	(u >   COUENNE_INFINITY) && (y0 > pow0)){ 

      if ((x0 > 0) && (y0 > pow0) ||  
	  (x0 < 0) && (y0 < pow0)) { // in orthant

	*brpts = 0;
	return fabs (pow0 - y0);

      } else {

	*brpts = x0;

	return (x0 > 0) ? // approx distance
	  projectSeg (x0, y0, x0, pow0, CoinMax (0., pow (y0, 1./k)), y0, 0) :
	  projectSeg (x0, y0, x0, pow0, CoinMin (0., pow (y0, 1./k)), y0, 0);
      }
    }

    // otherwise, the convexification is bounded.
    //
    // Apply min-area convexification

    if (u > l + COUENNE_EPS) {
      *brpts = safe_pow ((safe_pow (u, k) - safe_pow (l,k)) / (k* (u-l)), 1./(k-1.));
      if (u<0) 
	*brpts = -*brpts;
    }
  }

  if (k>1) { // case 5: k>1 /////////////////////////////////////////////////////////////////

    if (y0 < pow0) { // on the good side, i.e. out of the convex
		     // side. We don't care if u is infinity

      powertriplet pt (k);
      *brpts = powNewton (x0, y0, &pt);

      way = TWO_LEFT;	

      x0 -= *brpts;
      y0 -= pow (*brpts, k);

      return sqrt (x0*x0 + y0*y0);

    } else { // on the bad (concave) side, above curve

      // as a rule of thumb, take the x coordinate of the midpoint of
      // horizontal segment between current point and curve

      if (obj -> Strategy () == CouenneObject::MID_INTERVAL)
	*brpts = 0.5 * (x0 + pow (x0, 1. / k));
      else {
	powertriplet pt (k);
	*brpts = obj -> getBrPoint (&pt, x0, y0, l, u);
      }
      
      way = TWO_LEFT;

      if (l < 0.) l = 0.;

      CouNumber 
	powbpt = pow (*brpts, k),
	projL  = projectSeg (x0, y0, l, pow (l, k), *brpts, powbpt, -1);

      return (u > COUENNE_INFINITY) ? 
	CoinMin (projL, *brpts - x0) : 
	CoinMin (projL, projectSeg (x0, y0, *brpts, powbpt, u, pow (u,k), -1));
    }

  } else { // case 6: 0<k<1 //////////////////////////////////////////////////////////////////

    if (y0 < pow0) { // on the bad side, below

      // same rule of thumb as above, take the x coordinate of the
      // midpoint of horizontal segment between current point and
      // curve

      if (obj -> Strategy () == CouenneObject::MID_INTERVAL)
	*brpts = 0.5 * (x0 + pow (x0, 1. / k));
      else {
	powertriplet pt (k);
	*brpts = obj -> getBrPoint (&pt, x0, y0, l, u);
      }

      //*brpts = 0.5 * (x0 + pow (x0, 1. / k));
      way = TWO_LEFT;

      if (l < 0) l = 0;

      CouNumber 
	powbpt = pow (*brpts, k),
	projL  = projectSeg (x0, y0, l, pow (l, k), *brpts, powbpt, +1);

      return (u > COUENNE_INFINITY) ? 
	CoinMin (projL, powbpt - y0) : 
	CoinMin (projL, projectSeg (x0, y0, *brpts, powbpt, u, pow (u,k), +1));

    } else { // on the convex side. We don't care if u is infinity

      powertriplet pt (k);
      *brpts = powNewton (x0, y0, &pt);

      way = TWO_LEFT;

      x0 -= *brpts;
      y0 -= pow (*brpts, k);

      return sqrt (x0*x0 + y0*y0);
    }
  }

  // failsafe: return null index, so that CouenneObject::infeasibility()
  // picks the default variable/branchingpoint for the expression

  var = NULL;
  return 0.;
}
