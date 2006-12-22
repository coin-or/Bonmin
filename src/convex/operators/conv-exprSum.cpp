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


// Create standard formulation of sums of expressions

exprAux *exprSum::standardize (CouenneProblem *p) {

  // first of all, standardize all operands
  exprOp::standardize (p);

  // now simply return NULL, (the caller will assume there is nothing
  // to change), as a sum is already standard
  return p -> addAuxiliary (this);
}


// If this is a sum, no need to approximate from below or above but we
// rather add this as a linear constraints (all elements of arglist_
// are of the type c*x or x*c or x)

int exprSum::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  CouNumber sumconst = 0;

  nterms  = new int [1];
  *nterms = nargs_ + 1; // preliminary estimate of how many coefficient

  allocateCon (1, nterms, coeff, indices, rhs, sign);

  *nterms = 0;

  for (register int i=0; i<nargs_; i++)
    if (arglist_ [i] -> Type () == CONST)
      sumconst += arglist_ [i] -> Value ();
    else {
      convert_monomial (arglist_ [i], coeff [0] [*nterms], indices [0] [*nterms]);
      ++*nterms;
    }

  coeff   [0] [*nterms] = new exprConst (-1);
  indices [0] [*nterms] = w -> Index ();

  ++*nterms;

  *rhs  = new exprConst (- sumconst);
  *sign = COUENNE_EQ;

  return 1;
}


// generate convexification cut for constraint w = this

void exprSum::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

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

  printf ("Sum: "); cut -> print ();

  cs.insert (cut);
}
