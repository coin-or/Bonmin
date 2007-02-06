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


// generate convexification cut for constraint w = - x

void exprOpp::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {
  // easy... 

  if (cg -> isFirst ())
    cg -> addTangent (cs, w -> Index (), argument_ -> Index (), 0, 0, CouNumber (-1.), 0);
}
