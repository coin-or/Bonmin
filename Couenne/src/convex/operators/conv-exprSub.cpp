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

  int wi = w -> Index ();
  int xi = x -> Index ();
  int yi = y -> Index ();

  if (x->Type () == CONST) // (c - y) or (c - d)
    if (y->Type() == CONST) cg -> createCut (cs, x->Value()-y->Value(), 0, wi, 1, -1, 0, -1, 0, true);
    else                    cg -> createCut (cs, x->Value(),            0, wi, 1, yi, 1, -1, 0, true);
  else // (x - y) or (x - d)
    if (y->Type() == CONST) cg -> createCut (cs, y->Value(),     0, wi, -1., xi, 1., -1, 0.,true);
    else                    cg -> createCut (cs, 0.,          0, wi, -1., xi, 1., yi, -1.,true);
}
