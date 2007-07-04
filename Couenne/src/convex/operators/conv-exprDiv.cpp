/*
 * Name:    conv-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization and convexification methods for divisions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOp.hpp>
#include <exprDiv.hpp>
#include <exprClone.hpp>
#include <exprMul.hpp>
#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

// Create standard formulation of this expression
exprAux *exprDiv::standardize (CouenneProblem *p) {

  exprOp::standardize (p);
  return p -> addAuxiliary (this);
}


// generate convexification cut for constraint w = x/y

void exprDiv::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) {

  // TODO: Use method on Tawarmalani-Sahinidis

  // compute y bounds

  expression *yle, *yue, *ye = arglist_ [1];
  ye -> getBounds (yle, yue);

  CouNumber yl = (*yle) (), 
            yu = (*yue) (), k;

  delete yle;
  delete yue;

  int xi = arglist_ [0] -> Index (),
      wi = w  -> Index (),
      yi = ye -> Index ();

  bool cLW,  cRW,  cLY,  cRY = 
       cLW = cRW = cLY = true;

  if (!(cg -> isFirst ()) && chg) {
    cLW = chg [wi].lower != UNCHANGED;  cRW = chg [wi].upper != UNCHANGED;
    cLY = chg [yi].lower != UNCHANGED;  cRY = chg [yi].upper != UNCHANGED;
  }

  // if the denominator's bound interval has 0 as internal point,
  // there is no convexification

  if ((yl < -0) && (yu >  0)) return;

  // special case #1: y is almost constant (nonzero) --> y = k. We
  // only need a single plane w = x/k.

  if ((fabs (yl-yu) < COUENNE_EPS) && 
      ((k = fabs (yl+yu) / 2) > COUENNE_EPS)) {
    if (cLY || cRY) cg -> createCut (cs, 0., 0, wi, -1, xi, 1/k);
    return;
  }

  // compute w bounds

  expression *wle, *wue;
  w -> getBounds (wle, wue);

  CouNumber wl = (*wle) (), 
            wu = (*wue) ();

  delete wle; delete wue;

  // special case #2: w is almost constant (nonzero) --> w = x/y = k. We
  // only need a single plane x = y*k.

  if ((fabs (wl-wu) < COUENNE_EPS) &&
      ((k = fabs (wl+wu) / 2) > COUENNE_EPS)) {
    if (cLW || cRW) cg -> createCut (cs, 0., 0, yi, k, xi, -1.);
    return;
  }

  // Add McCormick convexification cuts. Reduce w = x/y to x = wy and
  // apply the same rule as for multiplications:
  //
  // 1) x >= yl w + wl y - yl wl
  // 2) x >= yu w + wu y - yu wu
  //
  // 3) x <= yl w + wu y - yl wu
  // 4) x <= yu w + wl y - yu wl

  if ((cLY || cLW) && (is_boundbox_regular (yl,wl))) cg->createCut (cs,yl*wl,-1,xi,-1.,wi,yl,yi,wl);
  if ((cRY || cRW) && (is_boundbox_regular (yu,wu))) cg->createCut (cs,yu*wu,-1,xi,-1.,wi,yu,yi,wu);
  if ((cLY || cRW) && (is_boundbox_regular (yl,wu))) cg->createCut (cs,yl*wu,+1,xi,-1.,wi,yl,yi,wu);
  if ((cRY || cLW) && (is_boundbox_regular (yu,wl))) cg->createCut (cs,yu*wl,+1,xi,-1.,wi,yu,yi,wl);
}
