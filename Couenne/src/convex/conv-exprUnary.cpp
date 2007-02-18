/*
 * Name:    conv-exprUnary.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to convexify n-ary operators
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprUnary.h>
#include <exprClone.h>

#include <CouenneProblem.h>


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)

exprAux *exprUnary::standardize (CouenneProblem *p) {

  exprAux *subst;

  if ((subst = argument_ -> standardize (p)))
    argument_ = new exprClone (subst);

  return p -> addAuxiliary (this);
}
