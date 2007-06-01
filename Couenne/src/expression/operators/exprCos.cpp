/*
 * Name:    exprCos.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for of cosine 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprCos.h>
#include <exprSin.h>
#include <exprBCos.h>
#include <exprOpp.h>
#include <exprMul.h>
#include <exprClone.h>

#include <math.h>


// return an expression -sin (argument), the derivative of cos (argument)
expression *exprCos::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprOpp (new exprSin (new exprClone (argument_)));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


// compute bounds of sin x given bounds of x 
void exprCos::getBounds (expression *&lb, expression *&ub) {
  //  lb = new exprConst (-1); 
  //  ub = new exprConst (1);
  //  return;

  expression *xl, *xu;
  argument_ -> getBounds (xl, xu);

  lb = new exprLBCos (xl, xu);
  ub = new exprUBCos (new exprClone (xl), new exprClone (xu));
}
