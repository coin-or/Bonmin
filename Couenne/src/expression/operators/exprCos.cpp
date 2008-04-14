/*
 * Name:    exprCos.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for cosines
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "exprCos.hpp"
#include "exprSin.hpp"
#include "exprBCos.hpp"
#include "exprOpp.hpp"
#include "exprMul.hpp"
#include "exprClone.hpp"


static const CouNumber 
  pi  = M_PI,
  pi2 = M_PI * 2.,
  pih = M_PI / 2.;


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


/// closest feasible points in function in both directions
void exprCos::closestFeasible (expression *varind, expression *vardep,
			       CouNumber& left, CouNumber& right) const
{
  CouNumber curr = (*varind)();
  int period = (int)(curr/pi2);
  CouNumber curr_noperiod = curr - pi2*period;
  CouNumber inv = acos((*vardep)());

  if (curr_noperiod < inv) {
    left = pi2*period - inv;
    right = pi2*period + inv;
  }
  else if (curr_noperiod < pi2-inv) {
    left = pi2*period + inv;
    right = pi2*(period+1) - inv;
  }
  else {
    left = pi2*(period+1) - inv;
    right = pi2*(period+1) + inv;
  }
}
