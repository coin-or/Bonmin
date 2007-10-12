/*
 * Name:    exprSin.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <exprSin.hpp>
#include <exprClone.hpp>
#include <exprCos.hpp>
#include <exprBSin.hpp>
#include <exprMul.hpp>


// differentiation

expression *exprSin::differentiate (int index) {

  return new exprMul (new exprCos (new exprClone (argument_)),
		      argument_ -> differentiate (index));
}


// compute bounds of sin x given bounds of x 

void exprSin::getBounds (expression *&lb, expression *&ub) {

  expression *xl, *xu;

  argument_ -> getBounds (xl, xu);

  lb = new exprLBSin (xl, xu);
  ub = new exprUBSin (new exprClone (xl), new exprClone (xu));
}
