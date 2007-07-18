/*
 * Name:    exprVar.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class for defining variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneCutGenerator.hpp>
#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprAux.hpp>
#include <exprOp.hpp>
#include <exprUnary.hpp>
#include <exprVar.hpp>
#include <exprBound.hpp>

#include <CouenneProblem.hpp>
#include <depGraph.hpp>


// is the variable one of those in varlist?

int exprVar::dependsOn (register int *varlist = NULL, register int n = 1) {

  if (!varlist) 
    return 1;

  while (n--) 
    if (varIndex_ == *varlist++) 
      return 1;

  return 0;
}


// Get lower and upper bound of a variable expression (if any)
void exprVar::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLowerBound (varIndex_); 
  ub = new exprUpperBound (varIndex_);
}


// generate convexification cut for constraint w = this
void exprVar::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) {
  if (cg -> isFirst ())
    cg -> createCut (cs, 0., 0, w -> Index (), 1., varIndex_);
}


/// implied bound processing. Expression w = x, upon change in lower
/// or upper bound of w, whose index is wind
bool exprVar::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  bool res;
  if (updateBound (-1, l + varIndex_, l [wind])) {res = true; chg [varIndex_].lower = CHANGED;}
  if (updateBound (+1, u + varIndex_, u [wind])) {res = true; chg [varIndex_].upper = CHANGED;}
  return res;
}


/// update dependence set with index of this variable
void exprVar::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) 
{dep -> insert (g -> lookup (varIndex_));}


/// Bound get
expression *exprVar::Lb () {return new exprLowerBound (varIndex_);}
expression *exprVar::Ub () {return new exprUpperBound (varIndex_);}
