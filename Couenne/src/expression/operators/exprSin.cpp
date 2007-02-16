/*
 * Name:    exprSin.C
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprSin.h>
#include <exprClone.h>
#include <exprCos.h>
#include <exprMul.h>


// differentiation

expression *exprSin::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprCos (new exprClone (argument_));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


// printing

void exprSin::print (std::ostream& out) const {
  exprUnary::print (out, "sin", PRE);
}
