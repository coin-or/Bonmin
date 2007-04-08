/*
 * Name:    exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for class logarithm
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <exprLog.h>
#include <exprConst.h>
#include <exprClone.h>
#include <exprMax.h>
#include <exprInv.h>
#include <exprMul.h>


/// get bounds of log (x) based on bounds of x

void exprLog::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  // [low|upp]er bound of w=log(x) is log (max (0, [low|upp]er (x)))
  lb = new exprLog (new exprMax (new exprConst (1e-100), lba));
  ub = new exprLog (new exprMax (new exprConst (1e-100), uba));
}


/// differentiation

expression *exprLog::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprInv (new exprClone (argument_));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


/// printing

void exprLog::print (std::ostream& out) const 
{exprUnary::print (out, "log", PRE);}


/// implied bound processing for expression w = log(x), upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprLog::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  int ind = argument_ -> Index ();

  bool res = updateBound (-1, l + ind, exp (l [wind]));
  res      = updateBound ( 1, u + ind, exp (u [wind])) || res;

  if (res)
    chg [ind] = 1;

  return res;
}
