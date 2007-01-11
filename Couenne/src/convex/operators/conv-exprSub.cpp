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

  // now simply return NULL, (the caller will assume there is nothing
  // to change), as subtraction is already standard
  return p -> addAuxiliary (this);
}


// if this is a subtraction, no need to approximate from below or
// above but we rather add this as a linear constraints (all elements
// of arglist_ are of the type c*x or x*c or x)
/*
int exprSub::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  nterms  = new int [1];
  *nterms = 3;

  allocateCon (1, nterms, coeff, indices, rhs, sign);

    // for each argument, two cases:
    //
    // 1) it is a multiplication c*x or x*c, hence we have to spot c
    //    and x and set c as coefficient and x as variable
    // 2) x is a variable, coefficient is one

  convert_monomial (arglist_ [0], coeff [0] [0], indices [0] [0]);
  convert_monomial (arglist_ [1], coeff [0] [1], indices [0] [1]);

  coeff [0] [1] = new exprOpp (coeff [0] [1]);

  coeff   [0] [2] = new exprConst (-1);
  indices [0] [2] = w -> Index ();

  *rhs  = new exprConst (0);
  *sign = COUENNE_EQ;

  return 1;
}
*/

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

  printf ("Sub: "); cut -> print ();

  cs.insert (cut);
}
