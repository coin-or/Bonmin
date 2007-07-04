/*
 * Name:    exprSin.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprSin.hpp>
#include <exprClone.hpp>
#include <exprCos.hpp>
#include <exprBSin.hpp>
#include <exprMul.hpp>


// differentiation

expression *exprSin::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprCos (new exprClone (argument_));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


// compute bounds of sin x given bounds of x 

void exprSin::getBounds (expression *&lb, expression *&ub) {

  //  lb = new exprConst (-1); 
  //  ub = new exprConst (1);
  //  return;

  expression *xl, *xu;

  argument_ -> getBounds (xl, xu);

  lb = new exprLBSin (xl, xu);
  ub = new exprUBSin (new exprClone (xl), new exprClone (xu));
}
