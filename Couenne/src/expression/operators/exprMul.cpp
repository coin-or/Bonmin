/*
 * Name:    exprMul.C
 * Author:  Pietro Belotti
 * Purpose: definition of multiplications
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprMul.h>
#include <exprSum.h>
#include <exprConst.h>
#include <CouennePrecisions.h>


// simplify multiplications

expression *exprMul::simplify () {

  exprOp:: simplify ();

  CouNumber prod = 1;

  bool found_one = false;

  for (register int i=0; i<nargs_; i++) {

    // check for null operands in multiplications

    if (arglist_ [i] -> Type () == CONST) {

      found_one = true;

      CouNumber c = arglist_ [i] -> Value ();
      prod *= c;

      if (fabs (c) < COUENNE_EPS_SIMPL) {

	for (int j=0; j<nargs_; j++)
	  if (arglist_ [j]) {
	    delete arglist_ [j];
	    arglist_ [j] = NULL;
	  }

	return new exprConst (0);
      }

      // check for nonzero constants in multiplications

      delete arglist_ [i];
      arglist_ [i] = NULL;
    }
  }

  if (found_one && shrink_arglist (prod, 1))
    return new exprConst (arglist_ [0] -> Value ());
  else return NULL;
}


// differentiate product of expressions

expression *exprMul::differentiate (int index) {

  expression **als = new expression * [nargs_];
  int nonconst = 0;

  for (int i = 0; i < nargs_; i++) 

    if (arglist_ [i] -> dependsOn (&index, 1)) {

      expression **alm = new expression * [nargs_];

      alm [i] = arglist_ [i] -> differentiate (index);

      for (int j = 0; j < nargs_; j++) 
	if (i!=j)
	  alm [j] = new exprCopy (arglist_ [j]);

      als [nonconst++] = new exprMul (alm, nargs_);
    }

  if (nonconst) 
    return new exprSum (als, nonconst);
  else {
    delete [] als;
    return new exprConst (0);
  }
}


// print

void exprMul::print (std::ostream& out)
  {exprOp::print (out, "*", INSIDE);}


// get a measure of "how linear" the expression is:
//
// CONSTANT  = 0: a constant
// LINEAR    = 1: linear
// QUADRATIC = 2: quadratic
// NONLINER  = 3: nonlinear non-quadratic

int exprMul::Linearity () {

  int lin = arglist_ [0] -> Linearity ();

  for (register int i=1; i<nargs_; i++) {

    lin += arglist_ [i] -> Linearity ();
    if (lin >= NONLINEAR) {
      lin = NONLINEAR;
      break;
    }
  }
  return lin;
}
