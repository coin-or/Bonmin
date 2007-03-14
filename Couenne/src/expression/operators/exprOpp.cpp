/*
 * Name:    exprOpp.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprOpp.h>


// find bounds of -x given bounds on x

void exprOpp::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  lb = new exprOpp (uba);
  ub = new exprOpp (lba);
}


// differentiation

inline expression *exprOpp::differentiate (int index) 
{return new exprOpp (argument_ -> differentiate (index));}


// printing

void exprOpp::print (std::ostream& out) const
{exprUnary::print (out, "-", PRE);}


/// implied bound processing for expression w = -x, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprOpp::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  int ind = argument_ -> Index ();

  bool res = updateBound (-1, l + ind, - u [wind]);
  res      = updateBound ( 1, u + ind, - l [wind]) || res;

  if (res)
    chg [ind] = 1;

  return res;
}
