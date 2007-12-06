/*
 * Name:    exprSin.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "exprSin.hpp"
#include "exprClone.hpp"
#include "exprCos.hpp"
#include "exprBSin.hpp"
#include "exprMul.hpp"


// differentiation

expression *exprSin::differentiate (int index) {

  return new exprMul (new exprCos (new exprClone (argument_)),
		      argument_ -> differentiate (index));
}


// compute bounds of sin x given bounds of x 

void exprSin::getBounds (expression *&lb, expression *&ub) {

  expression *xl, *xu;

  argument_ -> getBounds (xl, xu);

  lb = new exprLBSin (xl, xu);
  ub = new exprUBSin (new exprClone (xl), new exprClone (xu));
}


/// generalized implied bound procedure for sine/cosine
bool trigImpliedBound (enum cou_trig type, int wind, int xind,
		       CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  CouNumber *xl = l + xind, wl = l [wind],
            *xu = u + xind, wu = u [wind];

  bool tighter = false;

  /*
  //xl = new CouNumber;
  //xu = new CouNumber;
  for (int i = 0; i < 20; i++) {
  *xl = M_PI / 4. * (1. + i/5.);
  *xu = M_PI / 4. * (13. - i/3.);
  wl = -0.3;
  wu =  0.2;
  type = COU_COSINE;
  */

  CouNumber fl, fu, iwl, iwu, displacement;

  if (type == COU_SINE) {fl = sin (*xl); fu = sin (*xu); displacement = M_PI / 2.;} 
  else                  {fl = cos (*xl); fu = cos (*xu); displacement = 0;        }

  iwl = acos (wl);
  iwu = acos (wu);

  /*printf ("old bounds: [%g pi,%g pi] -> [%g,%g]  ---  w = [%g,%g] -8-> [%g pi, %g pi]\n", 
	  *xl / M_PI, *xu / M_PI, fl, fu, 
	  wl, wu, iwl / M_PI, iwu / M_PI);*/

  if (wu < fl) {

    CouNumber base = 2. * M_PI * floor ((*xl + M_PI/2) / (2.*M_PI)) + displacement;

    //printf ("wu, fl: base = %g pi\n", base / M_PI);

    if (updateBound (-1, xl, base + CoinMin (iwu, M_PI - iwu))) {
      tighter = true; 
      chg [xind]. setLower (t_chg_bounds::CHANGED);
    }
  }

  if (wu < fu) {

    CouNumber base = 2. * M_PI * floor ((*xu + M_PI/2)/ (2.*M_PI)) + displacement;

    //printf ("wu, fu: base = %g pi\n", base / M_PI);

    if (updateBound (+1, xu, base - CoinMin (iwu, M_PI - iwu))) {
      tighter = true; 
      chg [xind]. setUpper (t_chg_bounds::CHANGED);
    }
  }

  if (wl > fl) {

    CouNumber base = 2. * M_PI * floor ((*xl - M_PI/2) / (2.*M_PI)) + displacement + M_PI;

    //printf ("wl, fl: base = %g pi\n", base / M_PI);

    if (updateBound (-1, xl, base + CoinMin (iwl, M_PI - iwl))) {
      tighter = true; 
      chg [xind]. setLower (t_chg_bounds::CHANGED);
    }
  }

  if (wl > fu) {

    CouNumber base = 2. * M_PI * floor ((*xu - M_PI/2) / (2.*M_PI)) + displacement + M_PI;

    //printf ("wl, fu: base = %g pi\n", base / M_PI);

    if (updateBound (+1, xu, base - CoinMin (iwl, M_PI - iwl))) {
      tighter = true; 
      chg [xind]. setUpper (t_chg_bounds::CHANGED);
    }
  }

  //printf ("new bounds: [%g pi, %g pi]\n------------------------------\n", *xl / M_PI, *xu / M_PI);
  //}
  //exit (0);

  return tighter;
}
