/*
 * Name:    conv-exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of convexification methods for exprGroup
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprGroup.h>


/// Get lower and upper bound of an expression (if any)
void exprGroup::getBounds (expression *&lb, expression *&ub) {

}


/// reduce expression in standard form, creating additional aux
/// variables (and constraints)
exprAux *exprGroup::standardize (CouenneProblem *p) {

}


// generate equality between *this and *w
void exprGroup::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg) {

}
