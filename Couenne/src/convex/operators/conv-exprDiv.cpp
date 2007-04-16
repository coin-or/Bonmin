/*
 * Name:    conv-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization and convexification methods for divisions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOp.h>
#include <exprDiv.h>
#include <exprClone.h>
#include <exprMul.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

// Create standard formulation of this expression
exprAux *exprDiv::standardize (CouenneProblem *p) {

  exprOp::standardize (p);
  return p -> addAuxiliary (this);
}


// generate convexification cut for constraint w = x/y

void exprDiv::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // TODO: Use method on Tawarmalani-Sahinidis //////////////////////////////

  // get bounds of numerator and denominator

  expression *yle, *yue;

  arglist_ [1] -> getBounds (yle, yue);

  CouNumber yl = (*yle) (), 
            yu = (*yue) ();

  // if the denominator's bound interval has 0 as internal point,
  // there is no convexification

  if ((yl < - COUENNE_EPS) && 
      (yu >   COUENNE_EPS)) 
    return;

  expression *xle, *xue, *wle, *wue;

  arglist_ [0] -> getBounds (xle, xue);
  w            -> getBounds (wle, wue);

  expression *xe = arglist_ [0];
  expression *ye = arglist_ [1];

  CouNumber wl = (*wle) (), wu = (*wue) (),
            xl = (*xle) (), xu = (*xue) ();

  delete yle; delete yue;
  delete wle; delete wue;
  delete xle; delete xue;

  // Add McCormick convexification cuts. Reduce w = x/y to x = wy and
  // apply the same rule as for multiplications:
  //
  // 1) x >= yl w + wl y - yl wl
  // 2) x >= yu w + wu y - yu wu
  //
  // 3) x <= yl w + wu y - yl wu
  // 4) x <= yu w + wl y - yu wl

  int xi = xe -> Index (),
      wi = w  -> Index (),
      yi = ye -> Index ();

  if (is_boundbox_regular (yl, wl)) cg -> createCut (cs, yl*wl, -1, xi, -1., wi, yl, yi, wl);
  if (is_boundbox_regular (yu, wu)) cg -> createCut (cs, yu*wu, -1, xi, -1., wi, yu, yi, wu);
  if (is_boundbox_regular (yl, wu)) cg -> createCut (cs, yl*wu, +1, xi, -1., wi, yl, yi, wu);
  if (is_boundbox_regular (yu, wl)) cg -> createCut (cs, yu*wl, +1, xi, -1., wi, yu, yi, wl);
}
