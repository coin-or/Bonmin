/*
 * Name:    impliedBounds-exprDiv.cpp
 * Author:  Pietro Belotti
 * Purpose: implied bounds for division operators
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprDiv.h>

#include <CouennePrecisions.h>

/// positive part of a number x, [x]^+
inline CouNumber mypos (register CouNumber x)
{return ((x > 0.) ? x : 0.);}


/// negative part of a number x, [x]^-
inline CouNumber myneg (register CouNumber x)
{return ((x < 0.) ? x : 0.);}


/// implied bound processing for expression w = x/y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprDiv::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  bool res = false;

  // deal with the "y is a constant"
  if (arglist_ [1] -> Type () == CONST) {

    int ind = arglist_ [0] -> Index ();

    if (ind<0) {
      printf ("exprDiv::impliedBound: Warning, w=c/d constants\n");
      return false;
    }

    CouNumber c = arglist_ [1] -> Value ();

    if (fabs (c) < COUENNE_EPS) {
      printf ("exprDiv::impliedBound: Warning, division by zero\n");
      return false;
    }

    c = 1. / c;

    // a copy of exprMul::impliedBound for the case where y is a constant

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

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber x0;

    // deal with all other cases

    // Each bound on w is represented on the xy plane with two cones,
    // such that joining the extreme rays of both one obtains two
    // lines, or in other words, the second cone is obtained through
    // the transformation (x',y') = (-x,-y).
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
    // Where "INF" stands for "subproblem infeasible", "-" for
    // "nothing to improve", and the rest is improved (those with "?"
    // may improve).
 
    CouNumber *xl = l + xi, *yl = l + yi, wl = l [wind],
              *xu = u + xi, *yu = u + yi, wu = u [wind];

    // avoid changing bounds if x is constant
    if (xi == -1) 
      xl = xu = &x0;

    //////////// deal with lower bound of w=x/y /////////////////////////////////////////////

    if (wl < - COUENNE_EPS) { // columns 1 and 2 of 8 /////////////////////

      // (xl,yu) in upper and lower cone, or simply, yu > 0 and yu < 0. Row 2

      if ((*yu >   COUENNE_EPS) && (updateBound (+1, yu, mypos (*xl/wl)))) {res = true; chg [yi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (-1, xl, *yu * wl))) {res = true; chg [xi] = 1;}

      // (xu,yl), row 3

      if ((*yl >   COUENNE_EPS) && (updateBound (+1, xu, *yl * wl))) {res = true; chg [xi] = 1;}
      if ((*yl < - COUENNE_EPS) && (updateBound (-1, yl, myneg (*xu/wl)))) {res = true; chg [yi] = 1;}

    } else if (wl > COUENNE_EPS) { // columns 3 and 4 of 8 /////////////////

      // (xl,yl), row 1

      if ((*yl >   COUENNE_EPS) && (updateBound (-1, xl, *yl * wl))) {res = true; chg [xi] = 1;}
      if ((*yl < - COUENNE_EPS) && (updateBound (-1, yl, myneg (*xl/wl)))) {res = true; chg [yi] = 1;}

      // (xu,yu), row 4

      if ((*yu >   COUENNE_EPS) && (updateBound (+1, yu, mypos (*xu/wl)))) {res = true; chg [yi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (+1, xu, *yu * wl))) {res = true; chg [xi] = 1;}
    }

    ///////////// deal with upper bound of w=x/y ///////////////////////////////////////////

    if (wu < - COUENNE_EPS) { // columns 5 and 6 of 8 ///////////////////

      // (xl,yl), row 1

      if ((*yl >   COUENNE_EPS) && (updateBound (-1, xl, *yu * wu))) {res = true; chg [xi] = 1;}
      if ((*yl >   COUENNE_EPS) && (updateBound (-1, yl, *xu / wu))) {res = true; chg [yi] = 1;}
      if ((*yl < - COUENNE_EPS) && (*xl > *yl * wu) && 
	                           (updateBound (-1, yl, 0)))        {res = true; chg [yi] = 1;}

      // (xl,yu), row 2

      //    if ((*yu < - COUENNE_EPS) && (updateBound (+1, yu, *xl / wu))) {res = true; chg [yi] = 1;}

      // (xu,yl), row 3

      //    if ((*yl >   COUENNE_EPS) && (updateBound (-1, yl, *xu / wu))) {res = true; chg [yi] = 1;}

      // (xu,yu), row 4

      if ((*yu >   COUENNE_EPS) && (*xu < *yu * wu) && 
                                   (updateBound (+1, yu, 0)))        {res = true; chg [yi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (+1, xu, *yl * wu))) {res = true; chg [xi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (+1, yu, *xl / wu))) {res = true; chg [yi] = 1;}

    } else if (wu > COUENNE_EPS) { // columns 7 and 8 of 8 /////////////////////

      // (xl,yl), row 1

      //    if ((*yl >   COUENNE_EPS) && (updateBound (-1, yl, *xl / wu))) {res = true; chg [yi] = 1;}

      // (xl,yu), row 2

      if ((*yu >   COUENNE_EPS) && (*xl > *yu * wu) && 
                                   (updateBound (+1, yu, 0)))        {res = true; chg [yi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (+1, yu, *xu / wu))) {res = true; chg [yi] = 1;}
      if ((*yu < - COUENNE_EPS) && (updateBound (-1, xl, *yl * wu))) {res = true; chg [xi] = 1;}

      // (xu,yl), row 3

      if ((*yl >   COUENNE_EPS) && (updateBound (-1, yl, *xl / wu))) {res = true; chg [yi] = 1;}
      if ((*yl >   COUENNE_EPS) && (updateBound (+1, xu, *yu * wu))) {res = true; chg [xi] = 1;}
      if ((*yl < - COUENNE_EPS) && (*xu < *yl * wu) && 
                                   (updateBound (-1, yl, 0)))        {res = true; chg [yi] = 1;}

      // (xu,yu), row 4

      //    if ((*yu < - COUENNE_EPS) && (updateBound (?1, ??, *?? ? wu))) {res = true; chg [?i] = 1;}
    }
  }

  return res;
}
