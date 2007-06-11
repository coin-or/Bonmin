/*
 * Name:    exprPow.cpp
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


/// differentiate power of expressions

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


/// output

//void exprPow::print (std::ostream& out) const
//  {exprOp::print (out, "^", INSIDE);}


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


/// set implied bounds for function w = x^k, k negative, integer or
/// inverse integer, and odd

void invPowImplBounds (int, int, CouNumber *, CouNumber *, CouNumber, bool &, bool &);


/// implied bound processing for expression w = x^k, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprPow::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  bool resL, resU = resL = false;

  if (arglist_ [0] -> Type () <= CONST)   // base is constant or zero, nothing to do
    return false;

  if (arglist_ [1] -> Type () >  CONST) { // exponent is nonconstant, no good

    printf ("Couenne: Warning, power expression has nonconstant exponent: ");
    print (std::cout);
    printf ("\n");

    return false;
  }

  int index = arglist_ [0] -> Index ();

  CouNumber k = arglist_ [1] -> Value (); // exponent

  if ((fabs (k) < COUENNE_EPS) || 
      (fabs (k) > COUENNE_INFINITY)) // a null or infinite k is of little use
    return false;

  int intk; // integer (or integer inverse of) exponent

  bool isint    =           (k    - (intk = COUENNE_round (k))    < COUENNE_EPS), // k   is integer
       isinvint = !isint && (1./k - (intk = COUENNE_round (1./k)) < COUENNE_EPS); // 1/k is integer

  CouNumber wl = l [wind], // lower w
            wu = u [wind]; // upper w

  if (!(isint || isinvint) || (intk % 2)) { 

    // k or 1/k integer and odd, or non-integer --> it is a monotone
    // increasing function

    if (k > 0.) { // simple, just follow bounds

      /*resL = (wl > - COUENNE_INFINITY) ?
	updateBound (-1, l + index, pow (wl, 1./k)) : 
	updateBound (-1, l + index, - COUENNE_INFINITY);*/

      if (wl > - COUENNE_INFINITY) 
	resL = updateBound (-1, l + index, pow (wl, 1./k)); 


      /*resU = (wu < COUENNE_INFINITY) ?
	updateBound (+1, u + index, pow (wu, 1./k)) :
	updateBound (+1, u + index, COUENNE_INFINITY);*/

      if (wu < COUENNE_INFINITY)
	resU = updateBound (+1, u + index, pow (wu, 1./k));

    } else // slightly more complicated, resort to same method as in exprInv
      invPowImplBounds (wind, index, l, u, 1./k, resL, resU);
  } 
  else 
    if (isint) { // x^k, k integer and even --> non monotone

      CouNumber bound = (k<0) ? wl : wu;

      // |x| <= b^(1/k), where b is wl or wu depending on k negative
      // or positive, respectively

      if (bound > COUENNE_EPS) {

	if (fabs (bound) < COUENNE_INFINITY) {
	  resL = updateBound (-1, l + index, - pow (bound, 1./k));
	  resU = updateBound (+1, u + index,   pow (bound, 1./k));
	} /*else {
	  resL = updateBound (-1, l + index, - COUENNE_INFINITY);
	  resU = updateBound (+1, u + index,   COUENNE_INFINITY);
	  }*/
      }

      // invert check, if bounds on x do not contain 0 we may improve them

      bound = (k>0) ? wl : wu;

      CouNumber xl = l [index], 
	        xu = u [index],
                xb = pow (bound, 1./k);

      if      (xl > - xb + COUENNE_EPS) resL = updateBound (-1, l + index,   xb) || resL;
      else if (xu <   xb - COUENNE_EPS) resU = updateBound ( 1, u + index, - xb) || resU;

    } else { // x^k, k=(1/h), h integer and even, or x^k, neither k nor 1/k integer

      CouNumber lb = wl, ub = wu;

      if (k < 0) { // swap bounds as they swap on the curve x^k when 
	lb = wu;
	ub = wl;
      }

      if (lb > COUENNE_EPS) resL = updateBound (-1, l + index, pow (lb, 1./k));

      if (fabs (ub) < COUENNE_INFINITY) {
	if (ub > COUENNE_EPS) resU = updateBound (+1, u + index, pow (ub, 1./k));
      } //else                  resU = updateBound (+1, u + index, COUENNE_INFINITY);
    }

  if (resL) chg [index].lower = CHANGED;
  if (resU) chg [index].upper = CHANGED;

  return (resL || resU);
}
