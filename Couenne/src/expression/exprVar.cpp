/*
 * Name:    expression.C
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprBound.h>

#include <CouenneProblem.h>

// is the variable one of those in varlist?

bool exprVar::dependsOn (int *varlist = NULL, int n = 1) {

  if (!varlist) 
    return true;

  while (n--) 
    if (varIndex_ == *varlist++) 
      return true;

  return false;
}


// Get lower and upper bound of a variable expression (if any)

inline void exprVar::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLowerBound (varIndex_); 
  ub = new exprUpperBound (varIndex_);
}

// (possibly auxiliary) variable. Generate linear constraint w
// = x

int exprVar::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {
  nterms = new int;
  nterms [0] = 2;

  allocateCon (1, nterms, coeff, indices, rhs, sign);

  (*coeff) [0] = new exprConst (-1); 
  (*coeff) [1] = new exprConst (1);
  (*indices) [0] = w -> Index ();
  (*indices) [1] = varIndex_;
  *rhs = new exprConst (0);

  *sign = COUENNE_EQ;

  return 1;
}


// generate convexification cut for constraint w = this

void exprVar::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  OsiRowCut *cut = new OsiRowCut;

  CouNumber *coeff = new CouNumber [2];  
  int       *index = new int       [2];  

  cut -> setLb (0);
  cut -> setUb (0);

  coeff [0] =  1; index [0] = w -> Index ();
  coeff [1] = -1; index [1] = varIndex_;

  cs.insert (cut);
}
