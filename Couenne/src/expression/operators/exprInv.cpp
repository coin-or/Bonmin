/*
 * Name:    exprInv.C
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */


#include <exprInv.h>
#include <exprCopy.h>
#include <exprMul.h>


// differentiation

expression *exprInv::differentiate (int index) {

  expression **alm = new expression * [3];
  
  alm [0] = new exprInv (new exprCopy (argument_));
  alm [1] = new exprCopy (alm [0]);
  alm [2] = argument_ -> differentiate (index);

  return new exprMul (alm, 3);
}


// printing

void exprInv::print (std::ostream& out) {
  exprUnary::print (out, "1/", PRE);
}
