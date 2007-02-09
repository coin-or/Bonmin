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


// check if bounding box is suitable for a multiplication/division
// convexification constraint

bool is_boundbox_regular (CouNumber b1, CouNumber b2) {
  return 
    (fabs (b1) < COUENNE_INFINITY) && (fabs (b2) < COUENNE_INFINITY);
    //    && ((fabs (b1) > COUENNE_EPS) || (fabs (b2) > COUENNE_EPS));
}


// generate convexification cut for constraint w = x/y

void exprDiv::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

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

  CouNumber wl = (*wle) (), wu = (*wue) ();

  // Add McCormick convexification cuts. Reduce w = x/y to the case
  // x = wy and apply the same rule as for multiplications:
  //
  // 1) x >= yl w + wl y - yl wl
  // 2) x >= yu w + wu y - yu wu
  //
  // 3) x <= yl w + wu y - yl wu
  // 4) x <= yu w + wl y - yu wl

  int xi = xe -> Index (),
      wi = w  -> Index (),
      yi = ye -> Index ();

  OsiRowCut *cut;

  // set convexifier cuts global only if this is the broadest
  // convexification (that is, we are at the initial variable bounds).
  bool is_glob = cg -> isFirst ();

  // 1) 
  if (is_boundbox_regular (yl, wl)
      && (cut = cg -> createCut (yl*wl, -1, xi, CouNumber (-1.), wi, yl, yi, wl, is_glob)))
    cs.insert (cut);

  // 2) 
  if (is_boundbox_regular (yu, wu)
      && (cut = cg -> createCut (yu*wu, -1, xi, CouNumber (-1.), wi, yu, yi, wu, is_glob)))
    cs.insert (cut);

  // 3) 
  if (is_boundbox_regular (yl, wu)
      && (cut = cg -> createCut (yl*wu, +1, xi, CouNumber (-1.), wi, yl, yi, wu, is_glob)))
    cs.insert (cut);

  // 4) 
  if (is_boundbox_regular (yu, wl)
      && (cut = cg -> createCut (yu*wl, +1, xi, CouNumber (-1.), wi, yu, yi, wl, is_glob)))
    cs.insert (cut);
}
