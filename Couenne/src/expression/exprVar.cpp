/*
 * Name:    exprVar.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class for defining variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneCutGenerator.h>
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

void exprVar::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLowerBound (varIndex_); 
  ub = new exprUpperBound (varIndex_);
}


// generate convexification cut for constraint w = this

void exprVar::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {
  OsiRowCut *cut;

  if ((cut = cg -> createCut (0., 0, w -> Index (), 1., varIndex_)))
    cs.insert (cut);
}


/// implied bound processing. Expression w = x, upon change in lower
/// or upper bound of w, whose index is wind
bool exprVar::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  bool res = updateBound (-1, l + varIndex_, l [wind]);
  res      = updateBound (+1, u + varIndex_, u [wind]) || res;

  if (res) 
    chg [varIndex_] = 1;

  return res;
}
