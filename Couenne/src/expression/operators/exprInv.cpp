/*
 * Name:    exprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprInv.hpp"
#include "exprClone.hpp"
#include "exprMul.hpp"


// differentiation
expression *exprInv::differentiate (int index) {

  expression **alm = new expression * [3];
  
  alm [0] = new exprInv (new exprClone (argument_));
  alm [1] = new exprClone (alm [0]);
  alm [2] = argument_ -> differentiate (index);

  return new exprMul (alm, 3);
}


// printing
void exprInv::print (std::ostream &out, 
		     bool descend) const {
  out << "(1/";
  argument_ -> print (out, descend);
  out << ")";
}
//  exprUnary::print (out, "1/", PRE);}


/// general function to tighten implied bounds of a function w = x^k,
/// k negative, integer or inverse integer, and odd

void invPowImplBounds (int wind, int index, 
		       CouNumber *l, CouNumber *u, CouNumber k,
		       bool &resL, bool &resU) {

  CouNumber wl = l [wind], 
            wu = u [wind];

  // 0 <= l <= w <= u

  if (wl >= 0.) {
    if (wu > COUENNE_EPS) {
      if (wu < COUENNE_INFINITY) resL = updateBound (-1, l + index, pow (wu, k));
      else                       resL = updateBound (-1, l + index, 0.);
    }
    if (wl > COUENNE_EPS)        resU = updateBound (+1, u + index, pow (wl, k));
  }

  // l <= w <= u <= 0

  if (wu <= -0.) {
    if (wl < - COUENNE_EPS) {
      if (wl > - COUENNE_INFINITY) resU = updateBound (+1, u + index, pow (wl, k)) || resU;
      else                         resU = updateBound (+1, u + index, 0.)          || resU;
    }
    if (wu < - COUENNE_EPS)        resL = updateBound (-1, l + index, pow (wu, k)) || resL;
  }
}


/// implied bound processing for expression w = 1/x, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprInv::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  // Expression w = 1/x: we can only improve the bounds if 
  //
  //    0 <= l <= w <= u         or
  //         l <= w <= u <= 0. 
  //
  // Then 1/u <= x <= 1/l (given l, u finite and nonzero)

  int index = argument_ -> Index ();

  bool resL, resU = resL = false;

  invPowImplBounds (wind, index, l, u, -1., resL, resU);

  bool argInt = argument_ -> isInteger ();

  if (resL) {
    chg [index].setLower(t_chg_bounds::CHANGED);
    if (argInt) l [index] = ceil  (l [index] - COUENNE_EPS);
  }

  if (resU) {
    chg [index].setUpper(t_chg_bounds::CHANGED);
    if (argInt) u [index] = floor (u [index] + COUENNE_EPS);
  }

  return (resL || resU);
}
