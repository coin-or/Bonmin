/*
 * Name:    exprExp.C
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprExp.h>
#include <exprCopy.h>
#include <exprMul.h>


// differentiation

expression *exprExp::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprExp (new exprCopy (argument_));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


// printing

void exprExp::print (std::ostream& out)
  {exprUnary::print (out, "exp", PRE);}


// Get lower and upper bound of an expression (if any)
void exprExp::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  lb = new exprExp (lba);
  ub = new exprExp (uba);
}
