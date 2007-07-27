/*
 * Name:    exprCos.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for of cosine 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprCos.hpp>
#include <exprSin.hpp>
#include <exprBCos.hpp>
#include <exprOpp.hpp>
#include <exprMul.hpp>
#include <exprClone.hpp>

#include <math.h>


// return an expression -sin (argument), the derivative of cos (argument)
expression *exprCos::differentiate (int index) {

  return new exprOpp (new exprMul (new exprSin (new exprClone (argument_)),
				   argument_ -> differentiate (index)));
}


// compute bounds of sin x given bounds of x 
void exprCos::getBounds (expression *&lb, expression *&ub) {

  expression *xl, *xu;
  argument_ -> getBounds (xl, xu);

  lb = new exprLBCos (xl, xu);
  ub = new exprUBCos (new exprClone (xl), new exprClone (xu));
}
