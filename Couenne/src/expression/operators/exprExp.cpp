/*
 * Name:    exprExp.C
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprExp.h>
#include <exprClone.h>
#include <exprMul.h>


// differentiation

expression *exprExp::differentiate (int index) {

  expression **arglist = new expression * [2];

  arglist [0] = new exprExp (new exprClone (argument_));
  arglist [1] = argument_ -> differentiate (index);

  return new exprMul (arglist, 2);
}


// printing

void exprExp::print (std::ostream& out) const
  {exprUnary::print (out, "exp", PRE);}


// Get lower and upper bound of an expression (if any)
void exprExp::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  lb = new exprExp (lba);
  ub = new exprExp (uba);
}


/// implied bound processing for expression w = exp(x), upon change in
/// lower- and/or upper bound of w, whose index is wind
bool exprExp::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  bool res = false;
  int ind = argument_ -> Index ();

  CouNumber b;

  if ((b = l [wind]) >= COUENNE_EPS) // lower bound
    res = updateBound (-1, l + ind, log (b));

  if ((b = u [wind]) >= COUENNE_EPS) // upper bound
    res = updateBound ( 1, u + ind, log (b)) || res;

  if (res) 
    chg [ind] = 1;

  return res;
}
