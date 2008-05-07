/*
 * Name:    exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of divisions
 *
 * (C) Carnegie-Mellon University, 2006-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprDiv.hpp"
#include "exprConst.hpp"
#include "exprClone.hpp"
#include "exprMul.hpp"
#include "exprInv.hpp"
#include "exprSub.hpp"
#include "exprBDiv.hpp"

#include "CouennePrecisions.hpp"


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
	return new exprConst (0.);

      // otherwise, expression = k/y, return k*inv(y)

      expression *ret;

      if (fabs (arglist_ [0] -> Value ()-1) < COUENNE_EPS) {
	delete *arglist_;
	*arglist_ = NULL;
	ret = new exprInv (arglist_ [1]);
      }
      else ret = new exprMul (arglist_ [0], new exprInv (arglist_ [1]));

      arglist_ = NULL;
      return ret;
    }
  }
  else // only need to check if f2 == 0

    if (arglist_ [1] -> Type () == CONST) { // expression = x/h,
					    // transform into (1/h)*x

      expression *ret = new exprMul (new exprConst (1. / (arglist_ [1] -> Value ())), 
				     arglist_ [0]);
      delete arglist_ [1];
      arglist_ = NULL;
      return ret;
    }

  return NULL;
}


// differentiate quotient of expressions

expression *exprDiv::differentiate (int index) {

  if (!(arglist_ [0] -> dependsOn (index))  &&
      !(arglist_ [1] -> dependsOn (index)))
    return new exprConst (0.);

  expression **alm  = new expression * [2];
  expression **als  = new expression * [2];
  expression **alm2 = new expression * [3];

  exprInv *invg = new exprInv (new exprClone (arglist_ [1]));

  alm [0] = invg; // evaluated before alm2 [2]

  alm2 [0] = new exprClone (arglist_ [0]);
  alm2 [1] = arglist_ [1] -> differentiate (index);
  alm2 [2] = new exprClone (invg);

  als [0] = arglist_ [0] -> differentiate (index);
  als [1] = new exprMul (alm2, 3);

  alm [1] = new exprSub (als, 2);

  return new exprMul (alm, 2);
}


// get lower/upper bounds as a function of the arguments' lower/upper
// bounds

void exprDiv::getBounds (expression *&lb, expression *&ub) {

  expression **almin = new expression * [4];
  expression **almax = new expression * [4];

  arglist_ [0] -> getBounds (almin [0], almin [1]);
  arglist_ [1] -> getBounds (almin [2], almin [3]);

  almax [0] = new exprClone (almin [0]);
  almax [1] = new exprClone (almin [1]);
  almax [2] = new exprClone (almin [2]);
  almax [3] = new exprClone (almin [3]);

  lb = new exprLBDiv (almin, 4);
  ub = new exprUBDiv (almax, 4);
}


/// is this expression integer?
bool exprDiv::isInteger () {

  // only check if arguments (specifically, the denominator) are, *at
  // this point in the algorithm*, constant -- due to branching rules,
  // for instance. If so, check if the corresponding evaluated
  // expression is integer. Otherwise, check if denominator is +1 or
  // -1.

  expression *dl, *du, *nl, *nu;

  arglist_ [1] -> getBounds (dl, du);
  arglist_ [0] -> getBounds (nl, nu);

  register CouNumber 
    num = (*nl) (), 
    den = (*dl) ();

  bool
    denzero  = (fabs (den) < COUENNE_EPS),
    numconst = (fabs (num - (*nu) ()) < COUENNE_EPS);

  if ((fabs (num) < COUENNE_EPS)  && // numerator is zero
      numconst                    && // constant
      !denzero)                      // and denominator is nonzero

    return true;

  // otherwise...

  if (fabs (den - (*du) ()) < COUENNE_EPS) { // denominator is constant

    if (fabs (fabs (den) - 1) < COUENNE_EPS) // it is +1 or -1, check numerator
      return arglist_ [0] -> isInteger ();

    if (denzero) // it is zero, better leave...
      return false;

    if (numconst) { // numerator is constant, too

      CouNumber quot = num / den;

      if (fabs (COUENNE_round (quot) - quot) < COUENNE_EPS)
	return true;
    }
  }

  return false;
}


/// compute $y^{lv}$ and $y^{uv}$ for Violation Transfer algorithm
void exprDiv::closestFeasible (expression *varind,
			       expression *vardep, 
			       CouNumber &left,
			       CouNumber &right) const {

  expression *varoth = arglist_ [0]; // assume w = c/x

  bool numerator = false;

  if (varoth -> Index () == varind -> Index ()) { // actually w = x/c
    varoth = arglist_ [1];
    numerator = true;
  } else assert (arglist_ [1] -> Index () == varind -> Index ()); // right to assume w = c/x

  CouNumber 
    x = (*varind) (),
    y = (*vardep) (),
    c = (*varoth) ();

  if (numerator) // checking y = x/c

    if (c < 0.)
      if (y < x/c) right = c*y;
      else         left  = c*y;
    else 
      if (y < x/c) left  = c*y;
      else         right = c*y;

  else           // checking y = c/x

    if      (y < 0.) 
      if (x*y > c) right = c/y; // convex area in third orthant
      else         left  = c/y; // remaining of third+fourth orthant
    else if (y > 0.) 
      if (x*y > c) left  = c/y; // convex area in first orthant
      else         right = c/y; // remaining of first+second orthant
}


/// return l-2 norm of gradient at given point
CouNumber exprDiv::gradientNorm (const double *x) {

  int 
    ind0 = arglist_ [0] -> Index (),
    ind1 = arglist_ [1] -> Index ();

  CouNumber
    x0 = (ind0 < 0) ? fabs (arglist_ [0] -> Value ()) : fabs (x [ind0]),
    x1 = (ind1 < 0) ? fabs (arglist_ [1] -> Value ()) : fabs (x [ind1]),
    x1sq = x1 * x1;

  if (x1sq < 1/COUENNE_INFINITY) {
    x1sq = 1/COUENNE_INFINITY;
    if (x1 < 1/COUENNE_INFINITY) // implied
      x1 = 1/COUENNE_INFINITY;
  }

  if (ind0 < 0)
    if (ind1 < 0) return 0.;                // c/d
    else          return fabs (x0/(x1sq)); // c/y
  else 
    if (ind1 < 0) return 1. / x1;                                // x/d
    else          return sqrt (1. / x1sq + x0*x0 / (x1sq * x1sq)); // x/y
}

