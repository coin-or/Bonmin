/*
 * Name:    exprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of the absulute value of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprAbs.h>
#include <exprClone.h>
#include <exprMin.h>
#include <exprMax.h>
#include <exprOpp.h>


// find lower and upper bound of a given expression

void exprAbs::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;

  argument_ -> getBounds (lba, uba);

  expression **all = new expression * [6];
  expression **alu = new expression * [6];

  all [0] = new exprConst (0);  alu [0] = new exprConst (0);
  all [2] = new exprOpp (lba);  alu [2] = new exprOpp (new exprClone (lba));
  all [4] = uba;                alu [4] = new exprClone (uba);

  all [1] = new exprConst (0);  
  alu [1] = new exprMax (new exprOpp (new exprClone (lba)), new exprClone (uba));

  all [3] = new exprClone (lba); alu [3] = new exprClone (uba);
  all [5] = new exprOpp (new exprClone (uba));
  alu [5] = new exprOpp (new exprClone (lba));

  lb = new exprMin (all, 6);
  ub = new exprMin (alu, 6);
}


// differentiation

expression *exprAbs::differentiate (int index) {

  expression **arglist = new expression * [4];
  expression  *diffarg = argument_ -> differentiate (index);

  arglist [0] = new exprConst (0);
  arglist [1] = new exprClone (diffarg);
  arglist [2] = new exprOpp (new exprClone (argument_));
  arglist [3] = new exprOpp (diffarg);

  return new exprMin (arglist, 4);
}


// printing

void exprAbs::print (std::ostream& out) const {
  out << "|";
  argument_ -> print (out);
  out << "|";
}
