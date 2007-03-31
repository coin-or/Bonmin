/*
 * Name:    impliedBounds-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bounds for division operators
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprDiv.h>

#include <CouennePrecisions.h>


/// implied bound processing for expression w = x/y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprDiv::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  //  return false;

  bool resx, resy = resx = false;

  // deal with the "y is a constant"
  if (arglist_ [1] -> Type () == CONST) {

    int ind = arglist_ [0] -> Index ();

    if (ind < 0) {
      printf ("exprDiv::impliedBound: Warning, w=c/d constants\n");
      return false;
    }

    CouNumber c = arglist_ [1] -> Value ();

    if (fabs (c) < COUENNE_EPS) {
      printf ("exprDiv::impliedBound: Warning, division by zero\n");
      return false;
    }

    // a copy of exprMul::impliedBound for the case where y is a constant

    if (c > COUENNE_EPS) {

      resx = updateBound (-1, l + ind, l [wind] * c);
      resx = updateBound ( 1, u + ind, u [wind] * c) || resx;
    } 
    else if (c < - COUENNE_EPS) {

      resx = updateBound (-1, l + ind, u [wind] * c);
      resx = updateBound ( 1, u + ind, l [wind] * c) || resx;
    } 

    if (resx)
      chg [ind] = 1;

  } else {

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber x0;

    // deal with all other cases

    // Each bound on w is represented on the xy plane with two cones,
    // such that joining the extreme rays of both one obtains two
    // lines, or in other words, the second cone is obtained through
    // the transformation (x',y') = (-x,-y) applied to the first cone.
    //
    // Bounds can be tightened according to four different cases,
    // depending on which corner of the bounding box belongs to which
    // of the cones and on the sign of the bounds.
    //
    // Define wl <= w <= wu,
    //        xl <= x <= xu,
    //        yl <= y <= yu.
    //
    // Then the "tightenable" bounds are, depending on the corner:
    //
    //             _______________________________________________________
    //            |           w >= wl         |          w <= wu          |
    //            |___________________________|___________________________|
    //            |     l<0     |    l>0      |     u<0     |     u>0     |
    //            |_____________|_____________|_____________|_____________|
    //   Cone     |upper |lower |upper |lower |upper |lower |upper |lower |
    //            |      |      |      |      |      |      |      |      |
    // 1 xl,yl    | INF  |  -   |  xl  |  yl  |xl,yl?|  yl  |  yl  |  -   |
    // 2 xl,yu    |  yu  |  xl  |  -   | INF  |  -   |  yu  |  yu  |yu,xl?|
    // 3 xu,yl    |  xu  |  yl  | INF  |  -   |  yl  |  -   |yl,xu?|  yl  |
    // 4 xu,yu    |  -   | INF  |  yu  |  xu  |  yu  |xu,yu?|  -   |  yu  |
    //            |______|______|______|______|______|______|______|______|
    //
    // Where "INF" stands for "infeasible subproblem", "-" for
    // "nothing to improve", and the rest is improved (those with "?"
    // may improve).
 
    CouNumber *xl = l + xi, *yl = l + yi, wl = l [wind],
              *xu = u + xi, *yu = u + yi, wu = u [wind];

    /*printf ("from              : w[%d] [%e %e], x%d [%e %e] / y%d [%e %e]",
      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);*/

    // avoid changing bounds if x is constant
    if (xi == -1) 
      xl = xu = &x0;

    //////////// deal with lower bound of w=x/y /////////////////////////////////////////////

    if        (wl < - COUENNE_EPS) { // w >= wl, wl negative

      // point C: (xl,yl)

      resy = (*yl<0) && (*yl > *xl/wl) && updateBound (-1, yl, 0) || resy;//

      if ((*yl>0) && (*yl < *xl/wl)) { // point C violates x/y >= wl, down
	resx = updateBound (-1, xl, *yu*wl) || resx;//
	resy = updateBound (-1, yl, *xu/wl) || resy;//
      }

      // point B: (xu,yu)

      if ((*yu<0) && (*yu > *xu/wl)) { // point B violates x/y >= wl, down
	resx = updateBound (+1, xu, *yl*wl) || resx;
	resy = updateBound (+1, yu, *xl/wl) || resy;
      }

      resy = (*yu>0) && (*yu < *xu/wl) && updateBound (+1, yu, 0) || resy;

    } else if (wl >   COUENNE_EPS) { // w >= wl, wl positive

      //

      resy = (*yl<0) && (*yl < *xl/wl) && updateBound (-1, yl, mymin (*xl/wl, 0)) || resy;
      resx = (*yl>0) && (*yl > *xl/wl) && updateBound (-1, xl, *yl*wl)            || resx;

      //

      resx = (*yu<0) && (*yu < *xu/wl) && updateBound (+1, xu, *yu*wl)            || resx;
      resy = (*yu>0) && (*yu > *xu/wl) && updateBound (+1, yu, mymax (*xu/wl, 0)) || resy;
    }

    //////////// deal with upper bound of w=x/y /////////////////////////////////////////////

    if        (wu >   COUENNE_EPS) { // w <= wu, wu negative

      //

      resy = (*yl<0) && (*yl > *xu/wu) && updateBound (-1, yl, 0)      || resy;

      if ((*yl>0) && (*yl < *xu/wu)) {
	resx = updateBound (+1, xu, *yu*wu) || resx;
	resy = updateBound (-1, yl, *xl/wu) || resy;
      }

      //

      if ((*yu<0) && (*yu > *xl/wu)) {
	resx = updateBound (-1, xl, *yl*wu) || resx;
	resy = updateBound (+1, yu, *xu/wu) || resy;
      }

      resy = (*yu>0) && (*yu < *xl/wu) && updateBound (+1, yu, 0)      || resy;

    } else if (wu < - COUENNE_EPS) { // w <= wu, wu positive

      //

      resy = (*yl<0) && (*yl < *xu/wu) && updateBound (-1, yl, mymin (*xu/wu,0))  || resy;//
      resx = (*yu<0) && (*yu < *xl/wu) && updateBound (-1, xl, *yu*wu)            || resx;

      //

      resx = (*yl>0) && (*yl > *xu/wu) && updateBound (+1, xu, *yl*wu)            || resx;
      resy = (*yu>0) && (*yu > *xl/wu) && updateBound (+1, yu, mymax (*xl/wu,0))  || resy;
    }

    if (resx) chg [xi] = 1;
    if (resy) chg [yi] = 1;

    /*if (resx || resy) 
      printf ("                 \ntightened division: w[%d] [%e %e], x%d [%e %e] / y%d [%e %e]\n",
	      wind, wl, wu, xi, *xl, *xu, yi, *yl, *yu);
	      else printf ("                                                 \r");*/
  }

  return (resx || resy);
}
