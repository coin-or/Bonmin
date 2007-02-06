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


// generate convexification cut for constraint w = this

void exprSub::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  if (!(cg -> isFirst ()))
    return;

  CouNumber *coeff = new CouNumber [3];
  int       *index = new int       [3];
  OsiRowCut *cut   = new OsiRowCut;    

  CouNumber rhs = 0;
  int nvars = 1;

  coeff [0] = -1; 
  index [0] = w -> Index ();

  // first term

  if (arglist_ [0] -> Type () == CONST)
    rhs = - arglist_ [0] -> Value ();
  else {
    coeff [nvars] =  1; 
    index [nvars] = arglist_ [0] -> Index ();
    ++nvars;
  }

  // second term

  if (arglist_ [1] -> Type () == CONST)
    rhs += arglist_ [1] -> Value ();
  else {
    coeff [nvars] =  -1; 
    index [nvars] = arglist_ [1] -> Index ();
    ++nvars;
  }

  cut -> setUb (rhs);
  cut -> setLb (rhs);
  cut -> setRow (nvars, index, coeff);

  //  printf ("Sub: "); cut -> print ();

  cs.insert (cut);
}
