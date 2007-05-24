/*
 * Name:    conv-exprInv.cpp
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
  all [0] = new exprConst (0);       all [1] = new exprConst (- COUENNE_INFINITY);  // l<0<u
  all [2] = new exprOpp   (lba);     all [3] = new exprInv   (uba);                 // 0<l<u
  all [4] = new exprClone (uba);     all [5] = new exprInv   (new exprClone (uba)); // l<u<0

  lb = new exprMin (all, 6);

  expression **alu = new expression * [6];
  alu [0] = new exprConst (0);       alu [1] = new exprConst (COUENNE_INFINITY);   // l<0<u
  alu [2] = new exprClone (all [2]); alu [3] = new exprInv (new exprClone (lba));  // 0<l<u
  alu [4] = new exprClone (uba);     alu [5] = new exprInv (new exprClone (lba));  // l<u<0

  ub = new exprMin (alu, 6);
}


#define MIN_DENOMINATOR 1e-6

// generate convexification cut for constraint w = 1/x

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

  CouNumber x;

  int w_ind = aux       -> Index (), 
      x_ind = argument_ -> Index ();

  // special case: l and u are very close, replace function with
  // linear term

  if (fabs (u - l) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (u+l);
    cg -> createCut (cs, 2/x0, 0, w_ind, 1., x_ind, 1/(x0*x0));
    return;
  }

  // choose sampling points. If unbounded, bound using a rule of thumb

  int ns = cg -> nSamples ();

  if      (l < - COUENNE_INFINITY) l = ns * (u-1); // (-infinity, u] where u < 0
  else if (u >   COUENNE_INFINITY) u = ns * (l+1); // [l, +infinity) where l > 0

  // make bounds nonzero

  if (fabs (l) < COUENNE_EPS) {
    /*l /= (COUENNE_EPS / MIN_DENOMINATOR);
      if (fabs (l) < COUENNE_EPS) */
    l = (l<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;
  }

  if (fabs (u) < COUENNE_EPS) {
    /*u /= (COUENNE_EPS / MIN_DENOMINATOR);
      if (fabs (u) < COUENNE_EPS) */
      u = (u<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;
  }

  // bound
  cg -> addEnvelope 
    (cs, (l > 0) ? +1 : -1, 
     inv, oppInvSqr, w_ind, x_ind, 
     (cg -> isFirst ()) ? // is this first call?
       // place it somewhere in the interval (we don't care)
       ((l > COUENNE_EPS) ? l : u) :
       // otherwise, replace it where it gives deepest cut
       powNewton ((*argument_) (), (*aux) (), inv, oppInvSqr, inv_dblprime),
     l, u);

  delete xle; 
  delete xue;
}
