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

  if ((arglist_ [ind=0] -> Type () <= CONST) || 
      (arglist_ [ind=1] -> Type () <= CONST)) {

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

    /*printf ("from             : w[%d] [%e %e], x%d [%e %e] * y%d [%e %e]",
      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);*/

    // w's lower bound 

    bool resx = false, resy = false;

    if (wl >= 0.) {

      // point B in central infeasible area

      if (*xu * *yu < wl) {
	resx = (*xu * *yl < wl) && updateBound (+1, xu, wl / *yl) || resx;
	resy = (*xl * *yu < wl) && updateBound (+1, yu, wl / *xl) || resy;
      }

      // point C in central infeasible area

      if (*xl * *yl < wl) {
	resx = (*xl * *yu < wl) && updateBound (-1, xl, wl / *yu) || resx;
	resy = (*xu * *yl < wl) && updateBound (-1, yl, wl / *xu) || resy;
      }
    } else {

      // the infeasible set is a hyperbola with two branches

      // upper left
      resx = (*xl * *yl < wl) && (*yl > 0.) && updateBound (-1, xl, wl / *yl) || resx; // point C
      resy = (*xu * *yu < wl) && (*yu > 0.) && updateBound (+1, yu, wl / *xu) || resy; // point B

      // lower right
      resy = (*xl * *yl < wl) && (*yl < 0.) && updateBound (-1, yl, wl / *xl) || resy; // point C
      resx = (*xu * *yu < wl) && (*yu < 0.) && updateBound (+1, xu, wl / *yu) || resx; // point B
    }

    // w's upper bound 

    if (wu >= 0.) {

      // the infeasible set is a hyperbola with two branches

      // upper right
      resx = (*xu * *yl > wu) && (*yl > 0.) && updateBound (+1, xu, wu / *yl) || resx; // point D
      resy = (*xl * *yu > wu) && (*yu > 0.) && updateBound (+1, yu, wu / *xl) || resy; // point A

      // lower left
      resx = (*xl * *yu > wu) && (*yu < 0.) && updateBound (-1, xl, wu / *yu) || resx; // point A
      resy = (*xu * *yl > wu) && (*yl < 0.) && updateBound (-1, yl, wu / *xu) || resy; // point D

    } else {

      // point D in central infeasible area

      if (*xu * *yl > wu) {
	resx = (*xu * *yu > wu) && updateBound (+1, xu, wu / *yu) || resx;
	resy = (*xl * *yl > wu) && updateBound (-1, yl, wu / *xl) || resy;
      }

      // point A in central infeasible area

      if (*xl * *yu > wu) {
	resx = (*xl * *yl > wu) && updateBound (-1, xl, wu / *yl) || resx;
	resy = (*xu * *yu > wu) && updateBound (+1, yu, wu / *xu) || resy;
      }
    }

    /*if (resx || resy) 
      printf ("                 \ntightened product: w[%d] [%e %e], x%d [%e %e] * y%d [%e %e]\n",
	      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);
	      else printf ("                                                 \r");*/

    if (resx) chg [xi] = 1;
    if (resy) chg [yi] = 1;

    res = resx || resy;
  }

  return res;
}
