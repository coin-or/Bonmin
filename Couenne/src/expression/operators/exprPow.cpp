/*
 * Name:    exprPow.C
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouennePrecisions.h>
#include <exprPow.h>
#include <exprSum.h>
#include <exprMul.h>
#include <exprDiv.h>
#include <exprLog.h>
#include <exprConst.h>



// simplify power f(x) ^ g(x)

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
	return new exprConst (0);
  }
  else // only need to check if g(x) == 0

    if (arglist_ [1] -> Type () == CONST) {
      if (fabs (arglist_ [1] -> Value ()) < COUENNE_EPS_SIMPL) // expr = x ^ 0 = 1
	return new exprConst (1);
      else if (fabs (arglist_ [1] -> Value () - 1) < COUENNE_EPS_SIMPL) { // expr = x ^ 1 = x
	delete arglist_ [1];
	expression *ret = arglist_ [0];
	arglist_ [0] = arglist_ [1] = NULL;
	return ret;
      }
    }
  return NULL;
}


// differentiate power of expressions

expression *exprPow::differentiate (int index) {

  if (!(arglist_ [0] -> dependsOn (&index, 1))  &&
      !(arglist_ [1] -> dependsOn (&index, 1)))
    return new exprConst (0);

  expression **alm  = new expression * [2];
  expression **alp  = new expression * [2];
  expression **als  = new expression * [2];
  expression **alm1 = new expression * [2];
  expression **alm2 = new expression * [2];
  expression **ald  = new expression * [2];

  alp [0] = new exprCopy (arglist_ [0]);
  alp [1] = new exprCopy (arglist_ [1]);

  alm [0] = new exprPow (alp, 2);

  alm1 [0] = arglist_ [1] -> differentiate (index);
  alm1 [1] = new exprLog (new exprCopy (arglist_ [0]));

  als [0] = new exprMul (alm1, 2);

  ald [0] = new exprCopy (arglist_ [1]);
  ald [1] = new exprCopy (arglist_ [0]);

  alm2 [0] = new exprDiv (ald, 2);
  alm2 [1] = arglist_ [0] -> differentiate (index);

  als [1] = new exprMul (alm2, 2);

  alm [1] = new exprSum (als, 2);

  return new exprMul (alm, 2);
}


// output

void exprPow::print (std::ostream& out)
  {exprOp::print (out, "^", INSIDE);}


// get a measure of "how linear" the expression is:
//
// CONSTANT  = 0: a constant
// LINEAR    = 1: linear
// QUADRATIC = 2: quadratic
// NONLINER  = 3: nonlinear non-quadratic

int exprPow::Linearity () {

  if (arglist_ [0] -> Type () == CONST) {

    if (arglist_ [1] -> Type () == CONST) return CONSTANT;
    else                                  return NONLINEAR;
  }
  else {

    double exponent = arglist_ [1] -> Value ();

    if (fabs (exponent - FELINE_round (exponent)) > COUENNE_EPS)
      return NONLINEAR;

    if (arglist_ [1] -> Type () == CONST) { 

      int expInt = (int) FELINE_round (exponent);

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
