/*
 * Name:    conv-exprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to convexify an expression x^k, k constant
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.h>
#include <rootQ.h>
#include <exprPow.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprClone.h>
#include <exprMul.h>
#include <exprSum.h>
#include <exprLog.h>
#include <CouennePrecisions.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


std::map <int, CouNumber> Qroot::Qmap;

// Create standard formulation of this expression

exprAux *exprPow::standardize (CouenneProblem *p) {

  exprOp::standardize (p);

  if (arglist_ [0] -> Type () == CONST) { // expression is a^y, reduce
					  // to exp (x * log a)
    CouNumber base = arglist_ [0] -> Value ();

    if (fabs (base - M_E) < COUENNE_EPS)
      return p -> addAuxiliary (new exprExp (new exprClone (arglist_ [1])));
    else
      return p -> addAuxiliary 
	(new exprExp (p -> addAuxiliary (new exprMul (new exprClone (arglist_ [1]), 
						      new exprConst (log (base))))));
  }
  else

    if (arglist_ [1] -> Type () != CONST) // expression is x^y, reduce
					  // to exp (y*log(x));
      return p -> addAuxiliary 
	(new exprExp (p -> addAuxiliary 
		      (new exprMul 
		       (new exprClone (arglist_ [1]), 
			p -> addAuxiliary (new exprLog (new exprClone (arglist_ [0])))))));
    else                                  // expression is x^k, return
					  // as it is
	return p -> addAuxiliary (this);
}


// generate convexification cut for constraint w = x^k

void exprPow::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  return;

  // after standardization, all such expressions are of the form x^k,
  // with k constant

  CouNumber k = arglist_ [1] -> Value ();

  // get bounds of base

  expression *xle, *xue, 
    *wle, *wue, 
    *xe = arglist_ [0];

  xe  -> getBounds (xle, xue);
  aux -> getBounds (wle, wue);

  int w_ind = aux -> Index ();
  int x_ind = xe  -> Index ();

  CouNumber x = (*xe)  (), l = (*xle) (), u = (*xue) (),
            w = (*aux) ();

  OsiRowCut *cut;

  // classify power

  int intk = 0;

  bool isInt    =            fabs (k    - (intk = (int) (round (k))))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (intk = (int) (round (1./k)))) < COUENNE_EPS);
   
  // two macro-cases: 

  if (   (intk % 2) 
      && (k >   COUENNE_EPS) 
      && (l < - COUENNE_EPS) 
      && (u >   COUENNE_EPS)) {

    // 1) k (or its inverse) is positive, integer, and odd, and 0 is
    //    an internal point of the interval [l,u].

    Qroot qmap;

    // this case is somewhat simpler than the second, although we have
    // to resort to numerical procedures to find the (unique) root of
    // a polynomial Q(x) (see Liberti and Pantelides, 2003).

    CouNumber q = qmap (intk);

    int sign;

    if (isInvInt) {
      q = safe_pow (q, k);
      sign = -1;
    }
    else sign = 1;

    if (u > q * l) {
      addPowEnvelope (cg, cs, w_ind, x_ind, x, k, q*l, u, sign);
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), q*l, safe_pow (q*l,k), sign);
    }
    else
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), u, safe_pow (u,k), sign);

    // check if upper part needs a concave envelope

    if (l < q * u) {
      addPowEnvelope (cg, cs, w_ind, x_ind, x, k, l, q*u, -sign);
      cg -> addSegment (cs, w_ind, x_ind, q*u, safe_pow (q*u,k), u, safe_pow (u,k), -sign);
    }
    else
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), u, safe_pow (u,k), -sign);
  }
  else {

    // 2) all other cases.

    // if k is real or inv(k) is even, then lift l to max (0,l) but if
    // also u is negative, there is no convexification -- this
    // function is only defined on non-negative numbers

    if (!isInt 
	&& (!isInvInt || !(intk % 2))
	&& (l < - COUENNE_EPS) 
	&& (u < (l=0)))        // CAUTION! l is updated here, if negative
      return;

    // if k is negative and 0 is an internal point of [l,u], no
    // convexification is possible -- add a segment
    // between two tails of the asymptotes.

    if ((k < COUENNE_EPS) && (l < - COUENNE_EPS) && (u > COUENNE_EPS)) {

      if (!(intk % 2))
	cg -> addSegment (cs, w_ind, arglist_ [0] -> Index (), 
		    l, safe_pow (l,k), u, safe_pow (u,k), +1);
      return;
    }

    // ok, now we are ready to play. Between l and u we have a
    // convex/concave function that needs to be enveloped. Standard
    // segment and tangent cuts can be applied.

    // first of all, add trivial cut w >= 0 if violated

    if (!(intk % 2) && 
	(cg -> isFirst () || (w < - COUENNE_EPS)) &&
	(cut = cg -> createCut (0, +1, w_ind, CouNumber (1.), -1, 0., -1, 0., true)))
      cs.insert (cut);

    // create envelope. Choose sign based on k

    int sign = +1;

    // invert sign if (k is odd negative and l<0) or (k is in [0,1])
    if ((   l < - COUENNE_EPS) 
	&& (k < - COUENNE_EPS) 
	&& (intk % 2)
	|| (fabs (k-0.5) < 0.5 - COUENNE_EPS)
	|| ((u < - COUENNE_EPS) && (intk % 2) && (k > COUENNE_EPS)))
      sign = -1;

    // upper envelope -- when k negative, add only if bounds are far from 0

    if ((  (k > COUENNE_EPS)
	|| (l > COUENNE_EPS)
	|| (u < - COUENNE_EPS)) &&
	(l > - COUENNE_INFINITY + 1) &&
	(u <   COUENNE_INFINITY - 1)) {
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), u, safe_pow (u,k), - sign);
    }

    // similarly, pay attention not to add infinite slopes

#define POWER_RANGE 1e2;

    if (k > COUENNE_EPS) {

      if (u >   COUENNE_INFINITY - 1) u = x + POWER_RANGE;
      if (l < - COUENNE_INFINITY + 1) l = x - POWER_RANGE;
    }
    else {

      if (fabs (l) < COUENNE_EPS) l =  1. / POWER_RANGE; // l --> 0+
      if (fabs (u) < COUENNE_EPS) u = -1. / POWER_RANGE; // u --> 0-
    }

    addPowEnvelope (cg, cs, w_ind, x_ind, x, k, l, u, sign);
  }
}
