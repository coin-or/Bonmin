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

    if (fabs (base - M_E) < COUENNE_EPS) // is base e = 2.7182818...
      return p -> addAuxiliary (new exprExp (new exprClone (arglist_ [1])));
    else // no? convert k^x to e^(x log (k))
      return p -> addAuxiliary 
	(new exprExp (new exprClone 
		      (p -> addAuxiliary (new exprMul (new exprClone (arglist_ [1]), 
						       new exprConst (log (base)))))));
  }
  else
    if (arglist_ [1] -> Type () != CONST) //  x^y, convert to exp (y*log(x));
      return p -> addAuxiliary 
	(new exprExp (new exprClone (p -> addAuxiliary 
		      (new exprMul 
		       (new exprClone (arglist_ [1]), 
			new exprClone (p -> addAuxiliary 
				       (new exprLog (new exprClone (arglist_ [0])))))))));

    else return p -> addAuxiliary (this);  // expression is x^k, return as it is
}


// generate convexification cut for constraint w = x^k

void exprPow::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) {

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

  bool cL = !chg || (chg [x_ind].lower != UNCHANGED) || cg -> isFirst (),
       cR = !chg || (chg [x_ind].upper != UNCHANGED) || cg -> isFirst ();

  CouNumber w = (*aux) (), 
            x = (*xe)  (),
            l = (*xle) (), 
            u = (*xue) ();

  delete xle; delete xue;
  delete wle; delete wue;

  // if xl and xu are too close, approximate it as a line: sum the
  // segment through the two extreme points (l,l^k) and (u,u^k), and
  // the tangent at the midpoint ((l+u)/2, ((l+u)/2)^k)

  if (fabs (u-l) < COUENNE_EPS) {

    CouNumber avg     = 0.5 * (l+u),
              avg_k_1 = safe_pow (avg, k-1),
              lk      = safe_pow (l,   k),
              uk      = safe_pow (u,   k);

    if (cL || cR) cg -> createCut (cs, u*lk - l*uk + avg * avg_k_1 * (1-k), 0,
				   w_ind, u - l + 1, 
				   x_ind, lk-uk - k * avg_k_1);
    return;
  }

  // classify power

  int intk = 0;

  if (k < - COUENNE_INFINITY) { // w=x^{-inf} means w=0
    if (cL || cR) cg -> createCut (cs, 0., 0, w_ind, 1.);
    return;
  }

  if (k > COUENNE_INFINITY) // w=x^{inf} means not much...
    return;

  if (fabs (k) < COUENNE_EPS) { // w = x^0 means w=1
    if (cL || cR) cg -> createCut (cs, 1., 0, w_ind, 1.);
    return;
  }

  bool isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  // two macro-cases: 

  if (   (isInt || isInvInt)
      && (intk % 2) 
      && (k >   COUENNE_EPS) 
      && (l < - COUENNE_EPS) 
      && (u >   COUENNE_EPS)) {

    // 1) k (or its inverse) is positive, integer, and odd, and 0 is
    //    an internal point of the interval [l,u].

    // this case is somewhat simpler than the second, although if the
    // bound interval includes the origin we have to resort to
    // numerical procedures to find the (unique) root of a polynomial
    // Q(x) (see Liberti and Pantelides, 2003).

    CouNumber q = 1;

    if ((l<0) && (u>0)) {

      Qroot qmap;
      q = qmap (intk);
    }

    int sign;

    if (isInvInt) {
      if (cg -> isFirst ()) {
	w = (l>0) ? 1 : (u<0) ? -1 : 0;
	x = 0;
      }
      q = safe_pow (q, k);
      sign = -1;
    }
    else {
      if (cg -> isFirst ()) {
	x = (l>0) ? l : (u<0) ? u : 0;
	w = 0;
      }
      sign = 1;
    }

    // don't want big coefficients -- check only when k>1
    CouNumber powThres = (k<=1) ? COUENNE_INFINITY : pow (COU_MAX_COEFF, 1./k);

    // lower envelope
    if (l > -powThres) {
      if (l>0) addPowEnvelope (cg, cs, w_ind, x_ind, x, w, k,   l, u, sign); // 0<l<u, tangents only
      else if (u > q * l) { // upper x is after "turning point", add lower envelope
	addPowEnvelope        (cg, cs, w_ind, x_ind, x, w, k, q*l, u, sign);
	cg      -> addSegment     (cs, w_ind, x_ind, l, safe_pow (l,k), q*l, safe_pow (q*l,k), sign);
      } else cg -> addSegment     (cs, w_ind, x_ind, l, safe_pow (l,k), u,   safe_pow (u,  k), sign);
    }

    // upper envelope
    if (u < powThres) {
      if (u<0) addPowEnvelope (cg, cs, w_ind, x_ind, x, w, k, l,   u, -sign);  // l<u<0, tangents only
      else if (l < q * u) { // lower x is before "turning point", add upper envelope
	addPowEnvelope        (cg, cs, w_ind, x_ind, x, w, k, l, q*u, -sign);
	cg      -> addSegment     (cs, w_ind, x_ind, q*u, safe_pow (q*u,k), u, safe_pow (u,k), -sign);
      } else cg -> addSegment     (cs, w_ind, x_ind, l,   safe_pow (l,k),   u, safe_pow (u,k), -sign);
    }
  }
  else {

    // 2) all other cases.

    // if k is real or inv(k) is even, then lift l to max (0,l) but if
    // u is also negative, there is no convexification -- this
    // function is only defined on non-negative numbers

    if (!isInt 
	&& (!isInvInt || !(intk % 2))
	&& (l < - COUENNE_EPS) 
	&& (u < (l=0)))        // CAUTION! l updated if negative (k real)
      return;

    // if k is negative and 0 is an internal point of [l,u], no
    // convexification is possible -- just add a segment between two
    // tails of the asymptotes.

    if ((k < 0) && 
	(l < - COUENNE_EPS) && 
	(u >   COUENNE_EPS)) {

      if (!(intk % 2))
	cg -> addSegment (cs, w_ind, arglist_ [0] -> Index (), 
			  l, safe_pow (l,k), u, safe_pow (u,k), +1);
      return;
    }

    // Between l and u we have a convex/concave function that needs to
    // be enveloped. Standard segment and tangent cuts can be applied.

    // create upper envelope (segment)

    int sign = 1; // sign based on k

    // invert sign if 
    if (   ((l < - COUENNE_EPS) && (intk % 2) && (k < -COUENNE_EPS)) // k<0 odd, l<0
	|| ((u < - COUENNE_EPS) && (intk % 2) && (k >  COUENNE_EPS)) // k>0 odd, u<0
	|| (fabs (k-0.5) < 0.5 - COUENNE_EPS))                       // k in [0,1]
      sign = -1;

    CouNumber powThres = mymin (COUENNE_INFINITY, 
				pow (COU_MAX_COEFF, 1./k)), // don't want big coefficients
              powStep  = 1;

    // upper envelope
    if ((  (k > COUENNE_EPS)        // when k negative, add only if
	|| (l > COUENNE_EPS)        // bounds do not contain 0
	|| (u < - COUENNE_EPS)) &&
	(l > - powThres) &&         // and are finite
	(u <   powThres) &&
	(fabs (l+u) > COUENNE_EPS)) // bounds are not opposite (otherwise it's a variable bound)

      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l, k), u, safe_pow (u, k), -sign);

    // similarly, pay attention not to add infinite slopes

    if (cg -> isFirst()) {
      x = (k<0) ? ((u<0) ? u : (l>0) ? l : 0) : 0;
      w = 0;
    }

    if (k > COUENNE_EPS) {

      if (u >   powThres) u = x + powStep;
      if (l < - powThres) l = x - powStep;
    }
    else {

      if (fabs (l) < COUENNE_EPS) l =  1. / powThres; // l --> 0+
      if (fabs (u) < COUENNE_EPS) u = -1. / powThres; // u --> 0-
    }

    addPowEnvelope (cg, cs, w_ind, x_ind, x, w, k, l, u, sign);
  }
}
