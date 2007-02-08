/*
 * Name:    conv-exprSum.C
 * Author:  Pietro Belotti
 * Purpose: methods to standardize/convexify sum expressions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprSum.h>
#include <exprConst.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// Create standard formulation of sums of expressions

exprAux *exprSum::standardize (CouenneProblem *p) {

  // first of all, standardize all operands
  exprOp::standardize (p);

  // now simply return NULL, (the caller will assume there is nothing
  // to change), as a sum is already standard
  return p -> addAuxiliary (this);
}


// generate convexification cut for constraint w = this

void exprSum::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  if (!(cg -> isFirst ()))
    return;

  CouNumber *coeff = new CouNumber [nargs_ + 1];
  int       *index = new int       [nargs_ + 1];
  OsiRowCut *cut   = new OsiRowCut;

  CouNumber rhs = 0;

  coeff [0] = -1.; index [0] = w -> Index ();

  int nv = 1;

  for (int i=0; i<nargs_; i++) {

    if (arglist_ [i] -> Type () == CONST)
      rhs += arglist_ [i] -> Value ();
    else {
      coeff [nv] = 1.; 
      index [nv] = arglist_ [i] -> Index ();
      nv++;
    }
  }

  cut -> setRow (nv, index, coeff);
  cut -> setLb (-rhs);
  cut -> setUb (-rhs);

  delete [] index;
  delete [] coeff;

  cut -> setGloballyValid ();

  //  printf ("Sum: "); cut -> print ();

  cs.insert (cut);
}
