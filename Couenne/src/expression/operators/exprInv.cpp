/*
 * Name:    exprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */


#include <exprInv.h>
#include <exprClone.h>
#include <exprMul.h>


// differentiation

expression *exprInv::differentiate (int index) {

  expression **alm = new expression * [3];
  
  alm [0] = new exprInv (new exprClone (argument_));
  alm [1] = new exprClone (alm [0]);
  alm [2] = argument_ -> differentiate (index);

  return new exprMul (alm, 3);
}


// printing

void exprInv::print (std::ostream& out) const 
{exprUnary::print (out, "1/", PRE);}


/// general function (see below)
bool invPowImplBounds (int, int, CouNumber *, CouNumber *, CouNumber);


/// implied bound processing for expression w = 1/x, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprInv::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  // Expression w = 1/x: we can only improve the bounds if 
  //
  //    0 <= l <= w <= u         or
  //         l <= w <= u <= 0. 
  //
  // Then 1/u <= x <= 1/l (given l, u finite and nonzero)

  int index = argument_ -> Index ();

  bool res = invPowImplBounds (wind, index, l, u, -1.);

  if (res)
    chg [index] = 1;

  return res;
}


/// general function to tighten implied bounds of a function w = x^k,
/// k negative, integer or inverse integer, and odd

bool invPowImplBounds (int wind, int index, CouNumber *l, CouNumber *u, CouNumber k) {

  CouNumber wl = l [wind], 
            wu = u [wind];

  bool res = false;

  // 0 <= l <= w <= u

  if (wl >= 0.) {
    if (wu > COUENNE_EPS) {
      if (wu < COUENNE_INFINITY - 1) res = updateBound (-1, l + index, 1/wu);
      else                           res = updateBound (-1, l + index, 0.);
    }
    if (wl > COUENNE_EPS)            res = updateBound (+1, u + index, 1/wl) || res;
  }

  // l <= w <= u <= 0

  if (wu <= -0.) {
    if (wl < - COUENNE_EPS) {
      if (wl > - COUENNE_INFINITY + 1) res = updateBound (+1, u + index, 1/wl) || res;
      else                             res = updateBound (+1, u + index, 0.)   || res;
    }
    if (wu < - COUENNE_EPS)            res = updateBound (-1, l + index, 1/wu) || res;
  }

  return res;
}
