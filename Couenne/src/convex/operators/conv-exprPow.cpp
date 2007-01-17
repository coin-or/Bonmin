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
#include <exprDiv.h>
#include <exprOpp.h>
#include <exprSub.h>
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


// Find lower convex envelope for power x^k
/*
int exprPow::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // At this point, this should already be an expression of the sort
  // x^k, with x (auxiliary) variable and k constant. We need to deal
  // with ALL the following cases:
  //
  // 1) b   is integer and odd  (cube, x^5, etc)
  // 2) b   is integer and even (square, x^8, etc)
  // 3) 1/b is integer and odd  (cube root, x^(1/7), etc)
  // 4) 1/b is integer and even (square root, x^(1/4), etc)
  // 5) none of the above
  //
  // and with the case where the exponent is negative...

  if (arglist_ [1] -> Type () != CONST) {
    printf ("\n!!! ERROR !!! power is not in the form x^k\n\n");
    exit (-1);
  }

  CouNumber exponent = arglist_ [1] -> Value ();

  int rndexp;

  bool isInt    =  fabs (exponent - (rndexp = FELINE_round (exponent))) < COUENNE_EPS;
  bool isInvInt = !isInt &&  
    ((fabs (exponent) > COUENNE_EPS) && 
     (fabs (1/exponent - (rndexp = FELINE_round (1/exponent))) < COUENNE_EPS));

  // There are three possible convexifications:
  //
  // 1) segment + hull for the following cases:
  //
  //    a) k   positive, integer, even;
  //    b) 1/k integer, even;
  //    c) k   non-integer, both above and below 1;
  //
  // 2) segment + hull with checks for inclusion of 0 within bounds,
  // occurring when k<0 and k or 1/k integer and odd;
  //
  // 3) finally, when the expression is x^k or x^{1/k} with k integer
  // and odd, we have to add convexification depending on lower/upper
  // bound (see Liberti and Pantelides, 2003).
  //
  // Whether the segment and the linear hull are above or below
  // depends on the single case, but it suffices to change the signs.

  if (((exponent > 0) && isInt && !(rndexp % 2)) || // case 1a)
      (isInvInt && !(rndexp % 2)) ||                // case 1b)
      (!isInt && !isInvInt)) {                      // case 1c)

    nterms = new int [nSamples () + 1];
    for (int i=nSamples (); i>=0; i--)
      nterms [i] = 2;

    allocateCon (nSamples () + 1, nterms, coeff, indices, rhs, sign);

    expression **x_coeff = new expression * [nSamples ()];

    segment (coeff [0] [1], *rhs);
    hull    (x_coeff, rhs+1);

    enum con_sign hull_sign = COUENNE_GE,
                  seg_sign  = COUENNE_LE;

    // change signs if this is a concave rather than convex function
    // -- i.e., case b) and c) with k < 1

    if (isInvInt || (!isInt && (exponent < 1))) {
      hull_sign = COUENNE_LE;
      seg_sign  = COUENNE_GE;
    }

    int x_index = arglist_ [0] -> Index ();
    int w_index = w            -> Index ();

    // sets coefficients for hull part
    for (register int i = nSamples (); i; i--) {
      coeff   [i] [0] = new exprConst (1);        indices [i] [0] = w_index;
      coeff   [i] [1] = x_coeff [i-1];            indices [i] [1] = x_index;
      sign    [i]     = hull_sign;
    }

    coeff     [0] [0] = new exprConst (1);        indices [0] [0] = w_index;
    indices   [0] [1] = x_index;
    sign      [0]     = seg_sign;

    return nSamples () + 1;
  } else 
    if (exponent > 0) { // 3) integer odd, power or root 

    }
    else { // 2) need to check bounds

    }
    
  return 0;
}


// Method to create a linear, convex (concave) hull of a given
// function that is convex (concave)

void exprPow::hull (expression ** coeff, expression ** rhs) {

  expression *lba, *uba;
  arglist_ [0] -> getBounds (lba, uba);

  // now generate a set of underestimating constraints 

  int ns = nSamples ();

  ns--;

  for (int i=0; i<=ns; i++) {

    CouNumber alpha = (CouNumber) i / ns;

    expression *sample;

    if           (!i)  sample = uba;   // x_k = ub
    else if (i == ns)  sample = lba;   // x_k = lb
    else               sample =        // x_k = alpha lb + (1-alpha) ub, with 0 < alpha < 1

			 // introduce new operator?
      new exprSum
      (new exprMul    (new exprConst   (alpha), new exprClone (lba)),
       new exprMul    (new exprConst (1-alpha), new exprClone (uba)));

    expression *f_sample  = mirror   (sample);                 // value of f at sample: f  (x_k)
    expression *fp_sample = mirror_d (new exprClone (sample)); // derivative at x_k:    f' (x_k)

    coeff [i] = new exprOpp (new exprClone (fp_sample));

    rhs [i] = new exprSub (f_sample, new exprMul (fp_sample, new exprClone (sample)));
  }
}


// Method to create a segment approximating from above (below) a
// function that is convex (concave)

void exprPow::segment (expression*& coeff, expression*& rhs) {

  expression *lba, *uba;
  arglist_ [0] -> getBounds (lba, uba);

  // now generate one overestimating constraint

  // this is the slope of the segment over the curve exp (x) from lba
  // to uba: (exp(uba) - exp (lba)) / (uba - lba)
 
  expression *slope = new exprDiv (new exprSub (mirror        (lba), mirror        (uba)),
				   new exprSub (new exprClone (uba), new exprClone (lba)));
  coeff = slope;

  rhs = new exprDiv (new exprSub (new exprMul (mirror (lba), new exprClone (uba)),
				  new exprMul (mirror (uba), new exprClone (lba))),
		     new exprSub (new exprClone (uba), new exprClone (lba)));
}
*/

// generate convexification cut for constraint w = x^k

void exprPow::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // after standardization, all such expressions are of the form x^k,
  // with k constant

  CouNumber k = arglist_ [1] -> Value ();

  // get bounds of base

  expression *xle, *xue, 
             *wle, *wue;

  expression *xe = arglist_ [0];

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

    // check if lower part needs a convex envelope

    if (u > q * l) {
      addPowEnvelope (cg, cs, w_ind, x_ind, x, k, q*l, u, +1);
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), q, safe_pow (q,k), +1);
    }
    else
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), u, safe_pow (u,k), +1);

    // check if upper part needs a concave envelope

    if (l < q * u) {
      addPowEnvelope (cg, cs, w_ind, x_ind, x, k, l, q*u, -1);
      cg -> addSegment (cs, w_ind, x_ind, q, safe_pow (q,k), q, safe_pow (q,k), -1);
    }
    else
      cg -> addSegment (cs, w_ind, x_ind, l, safe_pow (l,k), u, safe_pow (u,k), -1);
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
	(cg -> isFirst () || 
	 !(cg -> addViolated ()) || 
	 (w < - COUENNE_EPS))) {

      cut = cg -> createCut (0, +1, w_ind, CouNumber (1.));

      //      if (cut) {printf ("Trivial cut: "); cut -> print ();}

      if (cut) 
	cs.insert (cut);
    }

    // create envelope. Choose sign based on k

    int sign = +1;

    // invert sign if (k is odd negative and l<0) or (k is in [0,1])
    if ((   l < - COUENNE_EPS) 
	&& (k < - COUENNE_EPS) 
	&& (intk % 2)
	|| (fabs (k-0.5) < 0.5 - COUENNE_EPS)
	|| ((u < - COUENNE_EPS) && (intk % 2) && (k > COUENNE_EPS)))
      sign = -1;

    // concave envelope -- when k negative, add only if bounds are far from 0

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
