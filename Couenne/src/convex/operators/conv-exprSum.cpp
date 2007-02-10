/*
 * Name:    conv-exprSum.cpp
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

  // TODO: Try a different approach. FLATTEN sum and reduce it to
  // something similar to a psg_elem: a sum of linear terms plus a
  // constant plus some nonlinear terms. (It only makes sense when
  // there is at least one monomial

  //  exprGroup *sum = new exprGroup (arglist_, nargs_);

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

  // first, make room for aux variable
  coeff [0] = -1.; index [0] = w -> Index ();

  int nv = 1;

  // scan arglist for (aux) variables and constants
  for (int i=0; i<nargs_; i++) {

    if (arglist_ [i] -> Type () == CONST)
      rhs += arglist_ [i] -> Value ();
    else {
      coeff [nv]   = 1.; 
      index [nv++] = arglist_ [i] -> Index ();
    }
  }

  cut -> setRow (nv, index, coeff);
  cut -> setUb (-rhs);
  cut -> setLb (-rhs);

  delete [] index;
  delete [] coeff;

  // added only once, it is global
  cut -> setGloballyValid ();

  cs.insert (cut);
}
