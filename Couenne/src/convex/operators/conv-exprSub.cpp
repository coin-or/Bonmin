/*
 * Name:    exprSub.C
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

  OsiRowCut *cut;

  expression *x = arglist_ [0];
  expression *y = arglist_ [1];

  int wind = w -> Index ();
  int xind = x -> Index ();
  int yind = y -> Index ();

  if (x->Type() == CONST) {

    CouNumber c = x -> Value ();

    if (y->Type() == CONST) cut = cg -> createCut (c - y->Value(), 0, wind, 1.);
    else                    cut = cg -> createCut (c,              0, wind, 1., yind, 1.);
  }
  else
    if (y->Type() == CONST) cut = cg -> createCut (y->Value(),     0, wind, -1., xind, 1.);
    else                    cut = cg -> createCut (0.,             0, wind, -1., xind, 1., yind, -1.);

  if (cut) {
    cut -> setGloballyValid ();
    cs.insert (cut);
  }
}
