/*
 * Name:    conv-exprInv.C
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the inverse operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>

#include <exprInv.h>
#include <exprClone.h>
#include <exprConst.h>
#include <exprMin.h>
#include <exprOpp.h>
#include <exprDiv.h>
#include <exprSum.h>
#include <exprMul.h>
#include <exprPow.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

// compute upper- and lower bounds of the expression w = 1/f(x) given
// bounds of f(x)

void exprInv::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  expression **all = new expression * [6];
  all [0] = new exprConst (0);                   all [1] = new exprConst (- COUENNE_INFINITY); 
  all [2] = new exprOpp   (new exprClone (lba)); all [3] = new exprInv   (new exprClone (uba)); 
  all [4] = new exprClone (uba);                 all [5] = new exprInv   (new exprClone (lba)); 

  lb = new exprMin (all, 6);

  expression **alu = new expression * [6];
  alu [0] = new exprConst (0);       alu [1] = new exprConst (COUENNE_INFINITY); 
  alu [2] = new exprClone (all [2]); alu [3] = new exprClone (all [5]); 
  alu [4] = new exprClone (uba);     alu [5] = new exprClone (all [3]); 

  ub = new exprMin (alu, 6);
}


// construct estimator for convex part of the inverse of an expression
/*
int exprInv::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {
  int ns = nSamples ();

  nterms = new int [ns];
  for (int i = ns; i--;)
    nterms [i] = 2;

  allocateCon (ns, nterms, coeff, indices, rhs, sign);

  expression *lb, *ub;
  argument_ -> getBounds (lb, ub);

  ns--;

  for (int i=0; i<=ns; i++) {

    expression *sample;

    if      (i==ns) sample = lb;
    else if (i==0)  sample = ub;
    else {

      CouNumber alpha = (CouNumber) i / ns;
      sample = new exprSum (new exprMul (new exprConst (alpha),   new exprClone (lb)),
			    new exprMul (new exprConst (1-alpha), new exprClone (ub)));
    }

    expression **alw = new expression * [6];
    expression **alx = new expression * [6];
    expression **alr = new expression * [6];

    // if min (0, -lb, ub) == 0 then bound interval contains 0 -> cannot
    // bound expression, return a null constraints

    alw [0] = new exprConst (0);   alw [1] = new exprConst (0); 
    alx [0] = new exprConst (0);   alx [1] = new exprConst (0);
    alr [0] = new exprConst (0);   alr [1] = new exprConst (0);

    // otherwise, if it is -lb, it means 0 <= lb <= ub and classical
    // convexification can be applied

    alw [2] = new exprOpp (new exprClone (lb));
    alx [2] = new exprOpp (new exprClone (lb));
    alr [2] = new exprOpp (new exprClone (lb));

    alw [3] = new exprConst (1);
    alx [3] = new exprPow (new exprClone (sample), new exprConst (-2));
    alr [3] = new exprDiv (new exprConst (2),      new exprClone (sample));

    // finally, if it is ub, then lb <= ub <= 0 and reverse
    // convexification is applied

    alw [4] = new exprClone (ub);
    alx [4] = new exprClone (ub);
    alr [4] = new exprClone (ub);

    alw [5] = new exprConst (-1);
    alx [5] = new exprOpp (new exprPow (new exprClone (sample), new exprConst (-2)));
    alr [5] = new exprDiv (new exprConst (-2), new exprClone (sample));

    // now specify coefficients of w and x and right-hand side

    coeff [i] [0] = new exprMin (alw, 6); indices [i] [0] = w -> Index ();
    coeff [i] [1] = new exprMin (alx, 6); indices [i] [1] = argument_ -> Index ();
    rhs   [i]     = new exprMin (alr, 6);
    sign  [i]     = COUENNE_GE;
  }

  return (ns+1);
}


// construct estimator for concave part of the inverse of an expression

int exprInv::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {
  nterms = new int [1];
  *nterms = 2;

  allocateCon (1, nterms, coeff, indices, rhs, sign);

  expression *lb, *ub;
  argument_ -> getBounds (lb, ub);

  expression **alw = new expression * [6];
  expression **alx = new expression * [6];
  expression **alr = new expression * [6];

  // if min (0, -lb, ub) == 0 then bound interval contains 0 -> cannot
  // bound expression, return a null constraints

  alw [0] = new exprConst (0);   alw [1] = new exprConst (0); 
  alx [0] = new exprConst (0);   alx [1] = new exprConst (0);
  alr [0] = new exprConst (0);   alr [1] = new exprConst (0);

  // otherwise, if it is -lb, it means 0 <= lb <= ub and classical
  // convexification can be applied

  alw [2] = new exprOpp (lb);
  alx [2] = new exprOpp (new exprClone (lb));
  alr [2] = new exprOpp (new exprClone (lb));

  alw [3] = new exprMul (new exprClone (lb), new exprClone (ub));
  alx [3] = new exprConst (1);
  alr [3] = new exprSum (new exprClone (lb), new exprClone (ub));

  // finally, if it is ub, then lb <= ub <= 0 and reverse
  // convexification is applied

  alw [4] = ub;
  alx [4] = new exprClone (ub);
  alr [4] = new exprClone (ub);

  alw [5] = new exprOpp (new exprMul (new exprClone (lb), new exprClone (ub)));
  alx [5] = new exprConst (1);
  alr [5] = new exprOpp (new exprSum (new exprClone (lb), new exprClone (ub)));

  // now specify coefficients of w and x and right-hand side

  coeff [0] [0] = new exprMin (alw, 6); indices [0] [0] = w -> Index ();
  coeff [0] [1] = new exprMin (alx, 6); indices [0] [1] = argument_ -> Index ();
  rhs   [0]     = new exprMin (alr, 6);
  sign  [0]     = COUENNE_LE;

  return 1;
}
*/

#define MIN_DENOMINATOR 1e-10

inline CouNumber oppInvSqr (CouNumber x) {
  CouNumber invx = inv (x); 
  return (- invx * invx);
}
//unary_function oppInvSqr;

// generate convexification cut for constraint w = this

void exprInv::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // get bounds of numerator and denominator

  expression *xle, *xue;

  argument_ -> getBounds (xle, xue);

  CouNumber l = (*xle) (), 
            u = (*xue) ();

  // if the denominator's bound interval has 0 as internal point,
  // there is no convexification

  if ((l < - COUENNE_EPS) && 
      (u >   COUENNE_EPS)) 
    return;

  expression *wle, *wue;

  aux -> getBounds (wle, wue);

  expression *xe = argument_;

  CouNumber x = (*xe)  (), xl = (*xle) (), xu = (*xue) (),
            w = (*aux) (), wl = (*wle) (), wu = (*wue) ();

  int w_ind = aux       -> Index (), 
      x_ind = argument_ -> Index ();

  // choose sampling
  // points. If unbounded, re-bound using a rule of thumb where

  int ns = cg -> nSamples ();

  if      (l < - COUENNE_INFINITY + 1) l = ns * u; // (-infinity, u] where u < 0
  else if (u >   COUENNE_INFINITY - 1) u = ns * l; // [l, +infinity) where l > 0

  // make bounds nonzero

  if (fabs (l) < COUENNE_EPS) {
    l /= (COUENNE_EPS / MIN_DENOMINATOR);
    if (fabs (l) < COUENNE_EPS) 
      l = (l<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;
  }

  if (fabs (u) < COUENNE_EPS) {
    u /= (COUENNE_EPS / MIN_DENOMINATOR);
    if (fabs (u) < COUENNE_EPS) 
      u = (u<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;
  }

  int sign = (l > 0) ? +1 : -1;

  // bound

  cg -> addEnvelope (cs, sign, inv, oppInvSqr, w_ind, x_ind, x, l, u);

  /*
  if ((cg -> ConvType () == UNIFORM_GRID) || cg -> isFirst ()) {

    // add tangent at each sampling point

    CouNumber sample = l, step = (u-l) / ns;

    for (int i = 0; i <= ns; i++) {

      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, inv (sample), - inv (sample*sample), sign);
      sample += step;
    }
  } else if (cg -> ConvType ()== CURRENT_ONLY)
    addTangent (cs, aux -> Index (), argument_ -> Index (), 
		x, inv (x), -inv(x*x), sign);
  else {

    CouNumber sample = x;

    addTangent (cs, aux -> Index (), argument_ -> Index (), 
		x, inv (x), -inv (x*x), sign);

    for (int i = 0; i <= ns/2; i++) {

      sample -= (x-l) / ns;
      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, inv (sample), - inv (sample*sample), sign);
    }

    sample = x;

    for (int i = 0; i <= ns/2; i++) {

      sample += (u-x) / ns;
      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, inv (sample), - inv (sample*sample), sign);
    }
  }
  */
}
