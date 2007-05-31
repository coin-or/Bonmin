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

bool exprMul::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  bool resL, resU = resL = false;
  int ind;

  if ((arglist_ [ind=0] -> Type () <= CONST) || 
      (arglist_ [ind=1] -> Type () <= CONST)) {

    // at least one constant in product w=cx:
    //
    // wl/c <= x <= wu/c, if c is positive
    // wu/c <= x <= wl/c, if c is negative

    CouNumber c = arglist_ [ind] -> Value ();

    // get the index of the nonconstant part
    ind = arglist_ [1-ind] -> Index ();

    if (ind==-1) // should not happen, it is a product of constant
      return false;

    if (c > COUENNE_EPS) {

      
      resL = (l [wind] > - COUENNE_INFINITY) && updateBound (-1, l + ind, l [wind] / c);
      resU = (u [wind] <   COUENNE_INFINITY) && updateBound ( 1, u + ind, u [wind] / c);
    } 
    else if (c < - COUENNE_EPS) {

      //      printf ("w_%d [%g,%g] = %g x_%d [%g,%g]\n", 
      //	      wind, l [wind], u [wind], c, ind, l [ind], u [ind]);

      resL = (u [wind] <   COUENNE_INFINITY) && updateBound (-1, l + ind, u [wind] / c);
      resU = (l [wind] > - COUENNE_INFINITY) && updateBound ( 1, u + ind, l [wind] / c);
    } 

    if (resL) chg [ind].lower = CHANGED;
    if (resU) chg [ind].upper = CHANGED;

      /*printf ("w_%d [%g,%g] -------> x_%d in [%g,%g] ", 
	      wind, l [wind], u [wind], 
	      ind,  l [ind],  u [ind]);*/
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

    bool resxL, resxU, resyL, resyU = 
      resxL = resxU = resyL = false;

    if (wl >= 0.) {

      // point B in central infeasible area

      if (*xu * *yu < wl) {
	resxU = (*xu * *yl < wl) && updateBound (+1, xu, wl / *yl);
	resyU = (*xl * *yu < wl) && updateBound (+1, yu, wl / *xl);
      }

      // point C in central infeasible area

      if (*xl * *yl < wl) {
	resxL = (*xl * *yu < wl) && updateBound (-1, xl, wl / *yu);
	resyL = (*xu * *yl < wl) && updateBound (-1, yl, wl / *xu);
      }
    } else {

      // the infeasible set is a hyperbola with two branches

      // upper left
      resxL = (*xl * *yl < wl) && (*yl > 0.) && updateBound (-1, xl, wl / *yl); // point C
      resyU = (*xu * *yu < wl) && (*yu > 0.) && updateBound (+1, yu, wl / *xu); // point B

      // lower right
      resyL = (*xl * *yl < wl) && (*yl < 0.) && updateBound (-1, yl, wl / *xl); // point C
      resxU = (*xu * *yu < wl) && (*yu < 0.) && updateBound (+1, xu, wl / *yu); // point B
    }

    // w's upper bound 

    if (wu >= 0.) {

      // the infeasible set is a hyperbola with two branches

      // upper right
      resxU = (*xu * *yl > wu) && (*yl > 0.) && updateBound (+1, xu, wu / *yl) || resxU; // point D
      resyU = (*xl * *yu > wu) && (*yu > 0.) && updateBound (+1, yu, wu / *xl) || resyU; // point A

      // lower left
      resxL = (*xl * *yu > wu) && (*yu < 0.) && updateBound (-1, xl, wu / *yu) || resxL; // point A
      resyL = (*xu * *yl > wu) && (*yl < 0.) && updateBound (-1, yl, wu / *xu) || resyL; // point D

    } else {

      // point D in central infeasible area

      if (*xu * *yl > wu) {
	resxU = (*xu * *yu > wu) && updateBound (+1, xu, wu / *yu) || resxU;
	resyL = (*xl * *yl > wu) && updateBound (-1, yl, wu / *xl) || resyL;
      }

      // point A in central infeasible area

      if (*xl * *yu > wu) {
	resxL = (*xl * *yl > wu) && updateBound (-1, xl, wu / *yl) || resxL;
	resyU = (*xu * *yu > wu) && updateBound (+1, yu, wu / *xu) || resyU;
      }
    }

    /*if (resx || resy) 
      printf ("                 \ntightened product: w[%d] [%e %e], x%d [%e %e] * y%d [%e %e]\n",
	      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);
	      else printf ("                                                 \r");*/

    if (resxL) chg [xi].lower = CHANGED;
    if (resxU) chg [xi].upper = CHANGED;
    if (resyL) chg [yi].lower = CHANGED;
    if (resyU) chg [yi].upper = CHANGED;

    resL = resxL || resyL;
    resU = resxU || resyU;
  }

  return (resL || resU);
}
