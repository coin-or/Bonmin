/*
 * Name:    conv-exprOp.C
 * Author:  Pietro Belotti
 * Purpose: methods to convexify n-ary operators
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOp.h>
#include <CouenneProblem.h>


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)
//
// For the base exprOp class we only do the first part (for argument
// list components only), and the calling class (Sum, Sub, Mul, Pow,
// and the like) will do the part for its own object

exprAux *exprOp::standardize (CouenneProblem *p) {

  exprVar *subst;

  for (register int i = nargs_; i--;)

    if ((subst = arglist_ [i] -> standardize (p)))
      arglist_ [i] = subst;

  return NULL;
}
