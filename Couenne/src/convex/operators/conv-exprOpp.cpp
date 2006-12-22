/*
 * Name:    conv-exprOpp.C
 * Author:  Pietro Belotti
 * Purpose: methods to convexify opposite of expressions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOpp.h>
#include <exprConst.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// construct linear under-estimator for expression within problem *p
// (p is used to add convexification constraints)

int exprOpp::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // we have an expression w = -x. We need to invert the linear hulls
  // and the sign of all appearances of

  nterms = new int [1];
  *nterms = 2;

  allocateCon (1, nterms, coeff, indices, rhs, sign);

  coeff [0] [0] = new exprConst (1); indices [0] [0] = w         -> Index ();
  coeff [0] [1] = new exprConst (1); indices [0] [1] = argument_ -> Index ();
  rhs   [0] = new exprConst (0);
  sign  [0] = COUENNE_EQ;

  return 1;
}


// generate convexification cut for constraint w = -x

void exprOpp::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // easy... 

  if (cg -> isFirst ())
    addTangent (cs, w -> Index (), argument_ -> Index (), 0, 0, -1., 0);
}
