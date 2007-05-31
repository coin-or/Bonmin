/*
 * Name:    exprSub.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of subtractions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprSub.h>
#include <exprOpp.h>
#include <CouennePrecisions.h>


// simplify subtractions

expression *exprSub::simplify () {

  exprOp:: simplify ();

  if ((*arglist_) -> Type () == CONST) { // expr = c1 - f2 

    CouNumber c0 = (*arglist_) -> Value ();

    if (arglist_ [1] -> Type () == CONST) { // expr = c1 - c2

      CouNumber c1 = arglist_ [1] -> Value ();

      delete arglist_ [0]; arglist_ [0] = NULL;
      delete arglist_ [1]; arglist_ [1] = NULL;

      return new exprConst (c0 - c1);
    }
    else if (fabs (c0) < COUENNE_EPS_SIMPL) { // expr = opp (f2)

      expression *ret = new exprOpp (arglist_ [1]);
      delete arglist_ [0];
      arglist_ [0] = arglist_ [1] = NULL;
      return ret;
    }
  }
  else // only need to check if f2 == 0

    if ((arglist_ [1] -> Type () == CONST) &&
	(fabs (arglist_ [1] -> Value ()) < COUENNE_EPS_SIMPL)) {
      // expr = f1 - 0 --> return f1

      expression *ret = arglist_ [0];
      delete arglist_ [1];
      arglist_ [0] = arglist_ [1] = NULL;
      return ret;
    }

  return NULL;
}


// differentiate product of expressions

expression *exprSub::differentiate (int index) {

  expression **arglist = new expression * [nargs_];

  for (int i = 0; i < nargs_; i++)
    if (arglist_ [i] -> dependsOn (&index, 1))
         arglist [i] = arglist_ [i] -> differentiate (index);
    else arglist [i] = new exprConst (0);

  return new exprSub (arglist, nargs_);
}


// print
//void exprSub::print (std::ostream& out) const
//  {exprOp::print (out, "-", INSIDE);}


// Get lower and upper bound of an expression (if any)
void exprSub::getBounds (expression *&lb, expression *&ub) {

  expression **alsl = new expression * [2];
  expression **alsu = new expression * [2];

  arglist_ [0] -> getBounds (alsl [0], alsu [0]);
  arglist_ [1] -> getBounds (alsu [1], alsl [1]);

  lb = new exprSub (alsl, 2);
  ub = new exprSub (alsu, 2);
}


/// implied bound processing for expression w = x-y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprSub::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  // caution, xi or yi might be -1
  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index ();

  if ((xi == -1) && (yi == -1)) // both x and y are constant
    return false;

  CouNumber xl, xu, yl, yu, 
            wl = l [wind], wu = u [wind];

  if (xi==-1) xl =         xu = arglist_ [0] -> Value ();
  else       {xl = l [xi]; xu = u [xi];}

  if (yi==-1) yl =         yu = arglist_ [1] -> Value ();
  else       {yl = l [yi]; yu = u [yi];}

  bool res = false;

  /// !!! REMOVED A return false;

  // w >= l

  if ((xi >= 0) && (updateBound (-1, l + xi, yl + wl))) {res = true; chg [xi].lower = CHANGED;}
  if ((yi >= 0) && (updateBound (+1, u + yi, xu - wl))) {res = true; chg [yi].upper = CHANGED;}

  // w <= u

  if ((xi >= 0) && (updateBound (+1, u + xi, yu + wu))) {res = true; chg [xi].upper = CHANGED;}
  if ((yi >= 0) && (updateBound (-1, l + yi, xl - wu))) {res = true; chg [yi].lower = CHANGED;}

  return res;
}
