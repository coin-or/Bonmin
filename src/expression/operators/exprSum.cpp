/*
 * Name:    exprSum.C
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprSum.h>
#include <exprCopy.h>
#include <exprConst.h>


// simplify sums

expression *exprSum::simplify () {

  exprOp:: simplify ();

  CouNumber total     = 0;
  bool   found_one = false;

  for (register int i=0; i<nargs_; i++) {

    // check for constant operands in multiplications

    if (arglist_ [i] -> Type () == CONST) {

      total += arglist_ [i] -> Value ();
      found_one = true;
      delete arglist_ [i]; 
      arglist_ [i] = NULL;
    }
  }

  if (found_one && shrink_arglist (total, 0))
    return new exprConst (arglist_ [0] -> Value ());
  else return NULL;
}


// differentiate sum of expressions

expression *exprSum:: differentiate (int index) {

  exprOp::differentiate (index);

  expression **arglist = new expression * [nargs_];

  register int nonconst = 0;

  for (int i = 0; i < nargs_; i++) 
    if (arglist_ [i] -> dependsOn (&index, 1))
      arglist [nonconst++] = arglist_ [i] -> differentiate (index);

  if (!nonconst) {
    delete [] arglist;
    return new exprConst (0);
  }
  else return new exprSum (arglist, nonconst);
}


// print

void exprSum::print (std::ostream& out)
  {exprOp::print (out, "+", INSIDE);}


// Get lower and upper bound of an expression (if any)

void exprSum::getBounds (expression *&lb, expression *&ub) {

  expression **all = new expression * [nargs_];
  expression **alu = new expression * [nargs_];

  for (register int i=0; i<nargs_; i++)
    arglist_ [i] -> getBounds (all [i], alu [i]);

  lb = new exprSum (all, nargs_);
  ub = new exprSum (alu, nargs_);
}


// get a measure of "how linear" the expression is:
//
// CONSTANT  = 0: a constant
// LINEAR    = 1: linear
// QUADRATIC = 2: quadratic
// NONLINER  = 3: nonlinear non-quadratic

int exprSum::Linearity () {

  int linmax = arglist_ [0] -> Linearity ();

  for (register int i=1; i<nargs_; i++) {
    register int lin = arglist_ [i] -> Linearity ();
    if (lin > linmax) linmax = lin;
  }
  return linmax;
}
