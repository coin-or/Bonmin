/*
 * Name:    exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for class logarithm
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "exprLog.hpp"
#include "exprConst.hpp"
#include "exprClone.hpp"
#include "exprMax.hpp"
#include "exprMin.hpp"
#include "exprDiv.hpp"


/// get bounds of log (x) based on bounds of x

void exprLog::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  // [low|upp]er bound of w=log(x) is log (max (0, [low|upp]er (x)))
  //  lb = new exprLog (new exprMax (new exprConst (1e-100), lba));
  //  ub = new exprLog (new exprMax (new exprConst (1e-100), uba));

  expression **all  = new expression * [4]; 

  all [0] = new exprClone (lba); all [1] = new exprLog (lba);
  all [2] = new exprConst (0);   all [3] = new exprConst (- COUENNE_INFINITY);
  lb = new exprMax (all, 4);

  expression **alu  = new expression * [4], 
             **alum = new expression * [4];

  alum [0] = new exprConst (COUENNE_INFINITY); 
  alum [1] = new exprConst (COUENNE_INFINITY); 
  alum [2] = new exprClone (uba); 
  alum [3] = new exprLog (uba);

  alu [0] = new exprClone (uba); alu [1] = new exprMin (alum, 4);
  alu [2] = new exprConst (0.);  alu [3] = new exprConst (- COUENNE_INFINITY);
  ub = new exprMax (alu, 4);
}


/// differentiation

expression *exprLog::differentiate (int index) {

  return new exprDiv (argument_ -> differentiate (index),
		      new exprClone (argument_));
}


/// implied bound processing for expression w = log(x), upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprLog::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  int ind = argument_ -> Index ();
  bool res = false;

  if (updateBound (-1, l + ind, exp (l [wind]))) {res = true; chg [ind].setLower(t_chg_bounds::CHANGED);}
  if (updateBound ( 1, u + ind, exp (u [wind]))) {res = true; chg [ind].setUpper(t_chg_bounds::CHANGED);}

  return res;
}
