/*
 * Name:    exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of divisions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprDiv.h>
#include <exprConst.h>
#include <exprMul.h>
#include <exprInv.h>
#include <exprSub.h>
#include <exprBDiv.h>

#include <CouennePrecisions.h>


// simplify division

expression *exprDiv::simplify () {

  exprOp:: simplify ();

  if ((*arglist_) -> Type () == CONST) { // expr = a / y 

    CouNumber c0 = (*arglist_) -> Value ();

    if (arglist_ [1] -> Type () == CONST) { // expr = a / b

      CouNumber c1 = arglist_ [1] -> Value ();

      delete arglist_ [0]; 
      delete arglist_ [1];
      arglist_ [0] = arglist_ [1] = NULL;

      return new exprConst (c0 / c1);
    }
    else {
      if (fabs (c0) < COUENNE_EPS_SIMPL) // expr = 0/y
	return new exprConst (0);

      // otherwise, expression = k/y, return k*inv(y)
      expression *ret = new exprMul (arglist_ [0], new exprInv (arglist_ [1]));
      arglist_ = NULL;
      return ret;
    }
  }
  else // only need to check if f2 == 0

    if (arglist_ [1] -> Type () == CONST) { // expression = x/h,
					    // transform into (1/h)*x

      expression *ret = new exprMul (new exprConst (1 / (arglist_ [1] -> Value ())), 
				     arglist_ [0]);
      delete arglist_ [1];
      arglist_ = NULL;
      return ret;
    }

  return NULL;
}


// differentiate quotient of expressions

expression *exprDiv::differentiate (int index) {

  if (!(arglist_ [0] -> dependsOn (&index, 1))  &&
      !(arglist_ [1] -> dependsOn (&index, 1)))
    return new exprConst (0);

  expression **alm  = new expression * [2];
  expression **als  = new expression * [2];
  expression **alm2 = new expression * [3];

  exprInv *invg = new exprInv (new exprCopy (arglist_ [1]));

  alm [0] = invg; // evaluated before alm2 [2]

  alm2 [0] = new exprCopy (arglist_ [0]);
  alm2 [1] = arglist_ [1] -> differentiate (index);
  alm2 [2] = new exprCopy (invg);

  als [0] = arglist_ [0] -> differentiate (index);
  als [1] = new exprMul (alm2, 3);

  alm [1] = new exprSub (als, 2);

  return new exprMul (alm, 2);
}


// output

void exprDiv::print (std::ostream& out) const
  {exprOp::print (out, "/", INSIDE);}


// get lower/upper bounds as a function of the arguments' lower/upper
// bounds

void exprDiv::getBounds (expression *&lb, expression *&ub) {

  expression **almin = new expression * [4];
  expression **almax = new expression * [4];

  arglist_ [0] -> getBounds (almin [0], almin [1]);
  arglist_ [1] -> getBounds (almin [2], almin [3]);

  almax [0] = new exprCopy (almin [0]);
  almax [1] = new exprCopy (almin [1]);
  almax [2] = new exprCopy (almin [2]);
  almax [3] = new exprCopy (almin [3]);

  lb = new exprLBDiv (almin, 4);
  ub = new exprUBDiv (almax, 4);
}
