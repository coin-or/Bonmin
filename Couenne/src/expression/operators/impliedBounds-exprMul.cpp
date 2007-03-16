/*
 * Name:    impliedBounds-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bounds for multiplications
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprMul.h>
#include <CouennePrecisions.h>


/// implied bound processing for expression w = x*y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprMul::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  bool res = false;
  int ind;

  if ((arglist_ [ind=0] -> Type () == CONST) || 
      (arglist_ [ind=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [ind] -> Value ();

    // get the index of the nonconstant part
    ind = arglist_ [1-ind] -> Index ();

    if (ind==-1) // should not happen, it is a product of constant
      return false;

    if (c > COUENNE_EPS) {

      res = updateBound (-1, l + ind, l [wind] / c);
      res = updateBound ( 1, u + ind, u [wind] / c) || res;

    } 
    else if (c < - COUENNE_EPS) {

      res = updateBound (-1, l + ind, u [wind] / c);
      res = updateBound ( 1, u + ind, l [wind] / c) || res;
    } 
    else res = false;

    if (res)
      chg [ind] = 1;

  } else {

    // these bounds would be implied by McCormick's convexification,
    // however we write them explicitly for internal use within bound
    // tightening, as otherwise they would only be known by Clp only.

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber *xl = l + xi, *yl = l + yi, wl = l [wind],
              *xu = u + xi, *yu = u + yi, wu = u [wind];

    // w's lower bound 

    if (wl > 0) {

      if ((*xl * *yl > wl) && 
	  (*xu * *yu < wl)) {

	if (updateBound (+1, xu, wl / *yl)) {res = true; chg [xi] = 1;}
	if (updateBound (+1, yu, wl / *xl)) {res = true; chg [yi] = 1;}
      } 
      else if ((*xl * *yl < wl) && 
	       (*xu * *yu > wl)) {

	if (updateBound (-1, xl, wl / *yu)) {res = true; chg [xi] = 1;}
	if (updateBound (-1, yl, wl / *xu)) {res = true; chg [yi] = 1;}
      }

    } else if (wl < 0) {

      if ((*xu * *yl > wl) && 
	  (*xl * *yu < wl)) {

	if (updateBound (-1, xl, wl / *yl)) {res = true; chg [xi] = 1;}
	if (updateBound (+1, yu, wl / *xu)) {res = true; chg [yi] = 1;}
      } 
      else if ((*xu * *yl < wl) && 
	       (*xl * *yu > wl)) {

	if (updateBound (+1, xu, wl / *yu)) {res = true; chg [xi] = 1;}
	if (updateBound (-1, yl, wl / *xl)) {res = true; chg [xi] = 1;}
      }
    }

    // w's upper bound 

    if (wu > 0) {

      if ((*xl * *yl > wu) && 
	  (*xu * *yu < wu)) {

	if (updateBound (-1, xl, wu / *yu)) {res = true; chg [xi] = 1;}
	if (updateBound (-1, yl, wu / *xu)) {res = true; chg [yi] = 1;}
      } 
      else if ((*xl * *yl < wu) &&
	       (*xu * *yu > wu)) {

	if (updateBound (+1, xu, wu / *yl)) {res = true; chg [xi] = 1;}
	if (updateBound (+1, yu, wu / *xl)) {res = true; chg [yi] = 1;}
      }

    } else if (wu < 0) {

      if ((*xu * *yl > wu) && 
	  (*xl * *yu < wu)) {

	if (updateBound (+1, xu, wu / *yu)) {res = true; chg [xi] = 1;}
	if (updateBound (+1, yl, wu / *xl)) {res = true; chg [yi] = 1;}
      } 
      else if ((*xu * *yl < wu) && 
	       (*xl * *yu > wu)) {

	if (updateBound (-1, xl, wu / *yl)) {res = true; chg [xi] = 1;}
	if (updateBound (+1, yu, wu / *xu)) {res = true; chg [yi] = 1;}
      }
    }
  }

  return res;
}
