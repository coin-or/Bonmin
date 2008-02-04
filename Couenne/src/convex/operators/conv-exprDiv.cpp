/*
 * Name:    conv-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization and convexification methods for divisions
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprAux.hpp"
#include "exprOp.hpp"
#include "exprDiv.hpp"
#include "exprClone.hpp"
#include "exprMul.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"


// Create standard formulation of this expression
exprAux *exprDiv::standardize (CouenneProblem *p, bool addAux) {

  exprOp::standardize (p);
  return (addAux ? (p -> addAuxiliary (this)) : new exprAux (this, p -> domain ()));
}


// generate convexification cut for constraint w = x/y
void exprDiv::generateCuts (expression *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  // compute y bounds

  CouNumber yl, yu;
  arglist_ [1] -> getBounds (yl, yu);

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index (),
      wi = w            -> Index ();

  bool cLW,  cRW,  cLY,  cRY = 
       cLW = cRW = cLY = true;

  if (!(cg -> isFirst ()) && chg) {
    cLW = chg [wi].lower() != t_chg_bounds::UNCHANGED;
    cRW = chg [wi].upper() != t_chg_bounds::UNCHANGED;
    cLY = chg [yi].lower() != t_chg_bounds::UNCHANGED;
    cRY = chg [yi].upper() != t_chg_bounds::UNCHANGED;
  }

  if ((yl < -0) && (yu > 0)) return;   // no convexification

  // special case #1: y is almost constant (nonzero) --> y = k. We
  // only need a single plane w = x/k.

  CouNumber k;

  if ((fabs (yl-yu) < COUENNE_EPS) && 
      ((fabs (k = yl+yu) / 2) > COUENNE_EPS)) {
    if (cLY || cRY)
      cg -> createCut (cs, 0., 0, wi, -1, xi, 1/k);
    return;
  }

  CouNumber wl, wu;
  w -> getBounds (wl, wu);

  if (lbw > wl) wl = lbw;
  if (ubw < wu) wu = ubw;

  // special case #2: w is almost constant (nonzero) --> w = x/y = k. We
  // only need a single plane x = y*k.

  if ((fabs (wl-wu) < COUENNE_EPS) &&
      ((k = fabs (wl+wu) / 2) > COUENNE_EPS)) {
    if (cLW || cRW)
      cg -> createCut (cs, 0., 0, yi, k, xi, -1.);
    return;
  }

  CouNumber xl, xu;
  arglist_ [0] -> getBounds (xl, xu);

  // same as product, just a change in coordinates

  CouNumber *x = w -> domain () -> x ();

  unifiedProdCuts (cg, cs,
		   wi, x [wi], wl, wu,
		   yi, x [yi], yl, yu,
		   xi, x [xi], xl, xu, chg);
}
