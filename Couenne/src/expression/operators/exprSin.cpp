/*
 * Name:    exprSin.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Carnegie-Mellon University, 2006-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "exprSin.hpp"
#include "exprClone.hpp"
#include "exprCos.hpp"
#include "exprBSin.hpp"
#include "exprMul.hpp"

static const CouNumber 
  pi  = M_PI,
  pi2 = M_PI * 2.,
  pih = M_PI / 2.;

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

  CouNumber fl, fu, iwl, iwu, displacement;

  if (type == COU_SINE) {fl = sin (*xl); fu = sin (*xu); displacement = pih;} 
  else                  {fl = cos (*xl); fu = cos (*xu); displacement = 0.;}

  iwl = acos (wl);
  iwu = acos (wu);

  /*printf ("### [%s] old bd: [%g pi,%g pi] -> [%g,%g]  ---  w = [%g,%g] -8-> [%g pi, %g pi]\n", 
	  type==COU_SINE ? "sin" : "cos", 
	   *xl / pi, *xu / pi, fl, fu, 
	   wl, wu, iwl/pi, iwu/pi);*/

  ////////////////////////////////////////////////////////////////////

  if (wu < fl - COUENNE_EPS) {

    CouNumber base = displacement + pi2 * floor ((*xl + pi - displacement) / pi2);

    //printf ("### wu, fl: base = %g pi\n", base / pi);

    if (updateBound (-1, xl, base + iwu)) {
      tighter = true; 
      chg [xind]. setLower (t_chg_bounds::CHANGED);
    }
  }

  if (wu < fu - COUENNE_EPS) {

    CouNumber base = displacement + pi2 * floor ((*xu + pi - displacement) / pi2);

    //printf ("### wu, fu: base = %g pi\n", base / pi);

    if (updateBound (+1, xu, base - iwu)) {
      tighter = true; 
      chg [xind]. setUpper (t_chg_bounds::CHANGED);
    }
  }

  if (wl > fl + COUENNE_EPS) {

    CouNumber base = displacement + pi2 * floor ((*xl - displacement) / pi2) + pi;

    //printf ("### wl, fl: base = %g pi\n", base / pi);

    if (updateBound (-1, xl, base + pi - iwl)) {
      tighter = true; 
      chg [xind]. setLower (t_chg_bounds::CHANGED);
    }
  }

  if (wl > fu + COUENNE_EPS) {

    CouNumber base = displacement + pi2 * floor ((*xu - displacement) / pi2) + pi;

    //printf ("### wl, fu: base = %g pi\n", base / pi);

    if (updateBound (+1, xu, base - pi + iwl)) {
      tighter = true; 
      chg [xind]. setUpper (t_chg_bounds::CHANGED);
    }
  }

  //printf ("### new bounds: [%g pi, %g pi]------------------------------\n", *xl/pi, *xu/pi);

  return tighter;
}

void exprSin::closestFeasible (expression *varind, expression *vardep,
			       CouNumber& left, CouNumber& right) const
{
  CouNumber curr = (*varind)() - pih;
  int period = (int)(curr/pi2);
  CouNumber curr_noperiod = curr - pi2*period;
  CouNumber inv = acos((*vardep)());

  if (curr_noperiod < inv) {
    left = pi2*period - inv;
    right = pi2*period + inv;
  }
  else if (curr_noperiod < pi2-inv) {
    left = pi2*period + inv;
    right = pi2*(period+1) - inv;
  }
  else {
    left = pi2*(period+1) - inv;
    right = pi2*(period+1) + inv;
  }
  left += pih;
  right += pih;
}
