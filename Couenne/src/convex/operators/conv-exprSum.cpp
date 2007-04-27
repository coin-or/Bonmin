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


/// generate convexification cut for constraint w = this

void exprSum::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  if (!(cg -> isFirst ()))
    return;

  CouNumber *coeff = new CouNumber [nargs_ + 1];
  int       *index = new int       [nargs_ + 1];
  OsiRowCut *cut   = new OsiRowCut;

  CouNumber rhs = 0;

  /// first, make room for aux variable
  coeff [0] = -1.; index [0] = w -> Index ();

  int nv = 1;

  /// scan arglist for (aux) variables and constants
  for (int i=0; i<nargs_; i++) {

    if (arglist_ [i] -> Type () == CONST)
      rhs -= arglist_ [i] -> Value ();
    else {
      coeff [nv]   = 1.; 
      index [nv++] = arglist_ [i] -> Index ();
    }
  }

  cut -> setRow (nv, index, coeff);
  cut -> setUb (rhs);
  cut -> setLb (rhs);

  delete [] index;
  delete [] coeff;

  /// added only once, it is global
  cut -> setGloballyValid ();

  cs.insert (cut);
}
