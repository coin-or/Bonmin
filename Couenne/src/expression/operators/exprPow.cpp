/*
 * Name:    exprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Carnegie-Mellon University, 2006-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>
#include <assert.h>

#include "CoinHelperFunctions.hpp"

#include "CouennePrecisions.hpp"
#include "exprPow.hpp"
#include "exprSum.hpp"
#include "exprMul.hpp"
#include "exprDiv.hpp"
#include "exprLog.hpp"
#include "exprConst.hpp"


/// simplify power f(x) ^ g(x)

expression *exprPow::simplify () {

  exprOp:: simplify ();

  if ((*arglist_) -> Type () == CONST) { // expr = c1 ^ g(x)

    CouNumber c0 = (*arglist_) -> Value ();

    if (arglist_ [1] -> Type () == CONST) { // expr = c1 ^ c2

      CouNumber c1 = arglist_ [1] -> Value ();

      delete arglist_ [0]; 
      delete arglist_ [1];

      arglist_ [0] = arglist_ [1] = NULL;

      return new exprConst (pow (c0, c1));
    }
    else 
      if (fabs (c0) < COUENNE_EPS_SIMPL) 
	return new exprConst (0.);
  }
  else // only need to check if g(x) == 0

    if (arglist_ [1] -> Type () == CONST) {

      CouNumber expon = arglist_ [1] -> Value ();

      if (fabs (expon) < COUENNE_EPS_SIMPL) // expr = x ^ 0 = 1
	return new exprConst (1.);

      else if (fabs (expon - 1) < COUENNE_EPS_SIMPL) { // expr = x ^ 1 = x

	delete arglist_ [1];
	expression *ret = arglist_ [0];
	arglist_ [0] = arglist_ [1] = NULL;
	return ret;
      }

      else if (fabs (expon + 1) < COUENNE_EPS_SIMPL) { // expr = x ^ -1 = 1/x

	delete arglist_ [1];
	expression *ret = new exprInv (arglist_ [0]);
	arglist_ [0] = arglist_ [1] = NULL;
	return ret;
      }
    }

  return NULL;
}


/// differentiate power of expressions

expression *exprPow::differentiate (int index) {

  if (!(arglist_ [0] -> dependsOn (&index, 1))  &&
      !(arglist_ [1] -> dependsOn (&index, 1)))
    return new exprConst (0.);

  // TODO: two cases

  expression **alm  = new expression * [2];
  expression **alp  = new expression * [2];
  expression **als  = new expression * [2];
  expression **alm1 = new expression * [2];
  expression **alm2 = new expression * [2];
  expression **ald  = new expression * [2];

  alp [0] = new exprClone (arglist_ [0]);
  alp [1] = new exprClone (arglist_ [1]);

  alm [0] = new exprPow (alp, 2);

  alm1 [0] = arglist_ [1] -> differentiate (index);
  alm1 [1] = new exprLog (new exprClone (arglist_ [0]));

  als [0] = new exprMul (alm1, 2);

  ald [0] = new exprClone (arglist_ [1]);
  ald [1] = new exprClone (arglist_ [0]);

  alm2 [0] = new exprDiv (ald, 2);
  alm2 [1] = arglist_ [0] -> differentiate (index);

  als [1] = new exprMul (alm2, 2);

  alm [1] = new exprSum (als, 2);

  return new exprMul (alm, 2);
}


/// get a measure of "how linear" the expression is:
///
/// ZERO      = 0: a zero
/// CONSTANT  = 1: a constant
/// LINEAR    = 2: linear
/// QUADRATIC = 3: quadratic
/// NONLINER  = 4: nonlinear non-quadratic

int exprPow::Linearity () {

  if (arglist_ [0] -> Type () == CONST) {

    if (arglist_ [1] -> Type () == CONST) return CONSTANT;
    else                                  return NONLINEAR;
  }
  else {

    double exponent = arglist_ [1] -> Value ();

    if (fabs (exponent - COUENNE_round (exponent)) > COUENNE_EPS)
      return NONLINEAR;

    if (arglist_ [1] -> Type () == CONST) { 

      int expInt = (int) COUENNE_round (exponent);

      if (arglist_ [0] -> Linearity () == LINEAR) {

	switch (expInt) {

	case 0:  return CONSTANT;
	case 1:  return LINEAR;
	case 2:  return QUADRATIC;

	default: return NONLINEAR;
	}
      }
      else 
	if (arglist_ [0] -> Linearity () == QUADRATIC) 
	  switch (expInt) {

	  case 0:  return CONSTANT;
	  case 1:  return QUADRATIC;

	  default: return NONLINEAR;
	  }
	else return NONLINEAR;
    }
    else return NONLINEAR;
  }
}


/// is this expression integer?
bool exprPow::isInteger () {

  // base

  if (!(arglist_ [0] -> isInteger ())) { 
    // base not integer: check if constant and integer
    CouNumber lb, ub;
    arglist_ [0] -> getBounds (lb, ub);

    if ((fabs (lb - ub) > COUENNE_EPS) ||
	!::isInteger (lb))
      return false;
  }

  // exponent

  if (!(arglist_ [1] -> isInteger ())) { 
    // exponent not integer: check if constant and integer (and
    // positive, or base negative integer)
    CouNumber lb, ub;
    arglist_ [1] -> getBounds (lb, ub);

    if ((fabs (lb - ub) > COUENNE_EPS) ||
	!::isInteger (lb))
      return false;

    if (lb < 0) { // exponent negative, check again base

      arglist_ [0] -> getBounds (lb, ub);

      if ((fabs (lb - ub) > COUENNE_EPS) ||
	  (fabs (lb) < COUENNE_EPS) ||
	  !::isInteger (1. / lb))
	return false;
    }
  }

  return true;
}


/// compute $y^{lv}$ and  $y^{uv}$ for Violation Transfer algorithm
void exprPow::closestFeasible (expression *varind,
			       expression *vardep, 
			       CouNumber &left,
			       CouNumber &right) const {
  CouNumber
    x  = (*varind) (), //info -> solution_ [indep     -> Index ()],
    y  = (*vardep) (), //info -> solution_ [dependent -> Index ()],
    k  = arglist_ [1] -> Value (),
    xk = safe_pow (x, k),
    yk = safe_pow (y, 1./k);

  int intk = 0;

  bool isInt    =            fabs (k    - (double) (intk = COUENNE_round (k)))    < COUENNE_EPS,
       isInvInt = !isInt && (fabs (1./k - (double) (intk = COUENNE_round (1./k))) < COUENNE_EPS);

  // three cases: 
  // 1) k or  1/k odd,        => have either left or right
  // 2) k or  1/k even,       => may have both
  // 3) k and 1/k fractional  => have either left or right

  if (isInt || isInvInt)

    if (intk % 2) // case 1

      if (k > 0) 
	((y < xk) ? left : right) = yk; // easy, x^k is continuous

      else

	if      (y < 0.)          // third, fourth orthant
	  if (y < xk) right = yk; // in convex region y < 1/x within third orthant
	  else        left  = yk; // remaining non-convex area

	else                      // first, second orthant
	  if (y > xk) left  = yk; // in convex region y > 1/x within first orthant
	  else        right = yk; // remaining non-convex area

    else // case 2

      if (y <= 0.) // third, fourth orthant => no solution
	left = - (right = COIN_DBL_MAX);

      else

	if (k > 0) 

	  if (k < 1) // roots, have x >= 0

	    if (x > yk) left  = yk;
	    else        right = yk;

	  else

	    if (x > yk)       left  =  yk;
	    else if (x < -yk) right = -yk;
	    else              left  = - (right = yk);

	else // k negative
	  if (y < xk) // between asymptotes
	    left = - (right = yk);
	  else  // in one of the two convex areas
	    if (x > 0) left  =  yk;
	    else       right = -yk;

  else // case 3: assume x bounded from below by 0

    if (k > 0) ((y < xk) ? left : right) = yk;
    else       ((y > xk) ? left : right) = yk;
}
