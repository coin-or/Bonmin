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


// convert an expression, which is either a multiplication of a
// constant and a variable or a variable, into a list of coefficient
// and indices to be fed into a LinearConstraint
/*
void convert_monomial (expression *arg, expression *&coeff, int &index) {

  // for each argument, three cases:
  //
  // 1) it is a multiplication c*x or x*c, hence we have to spot c
  //    and x and set c as coefficient and x as variable
  // 2) x is a variable, coefficient is one
  // 3) it is a constant
  if (arg -> Type () == N_ARY) {

    expression *term1 = arg -> ArgList () [0];
    expression *term2 = arg -> ArgList () [1];
      
    if (term1 -> Type () == CONST) {

      coeff = new exprConst (term1 -> Value ());
      index = term2 -> Index ();
    } 
    else {

      coeff = new exprConst (term2 -> Value ());
      index = term1 -> Index ();	
    }
  }
  else if (   (arg -> Type () == VAR) 
	   || (arg -> Type () == AUX)) {
    coeff = new exprConst (1);
    index = arg -> Index ();
  } else {
    printf ("CouError: found type %d instead of monomial\n", arg -> Type ());
    exit (-1);
  }
}
*/
