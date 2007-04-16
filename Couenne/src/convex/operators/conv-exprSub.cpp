/*
 * Name:    exprSub.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for the Subtraction class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <CouenneCutGenerator.h>
#include <exprSub.h>
#include <exprOpp.h>


// Create standard formulation of this expression

exprAux *exprSub::standardize (CouenneProblem *p) {

  // first of all, standardize all operands
  exprOp::standardize (p);

  return p -> addAuxiliary (this);
}


// generate convexification cut for constraint w = x - y

void exprSub::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  if (!(cg -> isFirst ()))
    return;

  // only add one cut at the beginning

  expression *x = arglist_ [0];
  expression *y = arglist_ [1];

  int wind = w -> Index ();
  int xind = x -> Index ();
  int yind = y -> Index ();

  if (x->Type () == CONST) { // (c - y) or (c - d)

    CouNumber c = x -> Value ();

    if (y->Type() == CONST) cg -> createCut (cs, c - y->Value(), 0, wind, 1.,   -1, 0., -1, 0., true);
    else                    cg -> createCut (cs, c,              0, wind, 1., yind, 1., -1, 0., true);
  }
  else // (x - y) or (x - d)
    if (y->Type() == CONST) cg -> createCut (cs, y->Value(),     0, wind, -1., xind, 1., -1, 0.,true);
    else                    cg -> createCut (cs, 0.,          0, wind, -1., xind, 1., yind, -1.,true);
}
