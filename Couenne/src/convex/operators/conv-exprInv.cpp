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


// derivative of inv (x)

inline CouNumber oppInvSqr (register CouNumber x) 
{return (- inv (x*x));}


#define MIN_DENOMINATOR 1e-10

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

  CouNumber x = (*xe)  ();

  int w_ind = aux       -> Index (), 
      x_ind = argument_ -> Index ();

  // choose sampling points. If unbounded, bound using a rule of thumb

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
}
