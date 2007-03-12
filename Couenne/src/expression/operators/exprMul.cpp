/*
 * Name:    exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of multiplications
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprMul.h>
#include <exprSum.h>
#include <exprConst.h>
#include <exprClone.h>
#include <CouennePrecisions.h>


/// simplify multiplications

expression *exprMul::simplify () {

  exprOp:: simplify ();

  if (nargs_ == 1) {

    expression *ret = arglist_ [0];
    arglist_ [0] = NULL;
    return ret;
  }

  CouNumber prod = 1;

  bool found_one = false;

  for (register int i=0; i<nargs_; i++) {

    // check for null operands in multiplications

    if (arglist_ [i] -> Type () == CONST) {

      found_one = true;

      CouNumber c = arglist_ [i] -> Value ();
      prod *= c;

      if (fabs (c) < COUENNE_EPS_SIMPL) {

	for (int j=0; j<nargs_; j++)
	  if (arglist_ [j]) {
	    delete arglist_ [j];
	    arglist_ [j] = NULL;
	  }

	return new exprConst (0);
      }

      // check for nonzero constants in multiplications

      delete arglist_ [i];
      arglist_ [i] = NULL;
    }
  }
  /*
  if (found_one && shrink_arglist (prod, 1))
    return new exprConst (arglist_ [0] -> Value ());
  */
  if (found_one && shrink_arglist (prod, 1)) {
    expression *ret = arglist_ [0];
    arglist_ [0] = NULL;
    return ret;
  }
  else return NULL;
}


// differentiate product of expressions

expression *exprMul::differentiate (int index) {

  expression **als = new expression * [nargs_];
  int nonconst = 0;

  for (int i = 0; i < nargs_; i++) 

    if (arglist_ [i] -> dependsOn (&index, 1)) {

      expression **alm = new expression * [nargs_];

      alm [i] = arglist_ [i] -> differentiate (index);

      for (int j = 0; j < nargs_; j++) 
	if (i!=j)
	  alm [j] = new exprClone (arglist_ [j]);

      als [nonconst++] = new exprMul (alm, nargs_);
    }

  if (nonconst) 
    return new exprSum (als, nonconst);
  else {
    delete [] als;
    return new exprConst (0);
  }
}


/// print

void exprMul::print (std::ostream& out) const
  {exprOp::print (out, "*", INSIDE);}


// get a measure of "how linear" the expression is:
//
// ZERO      = 0: constant 0
// CONSTANT  = 1: a constant
// LINEAR    = 2: linear
// QUADRATIC = 3: quadratic
// NONLINER  = 4: nonlinear non-quadratic

int exprMul::Linearity () {

  int lin0 = arglist_ [0] -> Linearity ();

  if (lin0 >= NONLINEAR) return NONLINEAR;
  if (lin0 == ZERO)      return ZERO;

  for (register int i=1; i<nargs_; i++) {

    int lin = arglist_ [i] -> Linearity ();

    switch (lin) {
    case NONLINEAR: return NONLINEAR;
    case ZERO:      return ZERO;
    case LINEAR:    lin0++; break;
    case QUADRATIC: lin0 += 2; break;
    default: break;
    }

    if (lin0 >= NONLINEAR) return NONLINEAR;
  }
  return lin0;
}


// choose which, between x and y, to branch on. Same choice as in
// exprMul, this function is defined in branch/getFixVarBinFun.cpp
expression *getFixVarBinFun (expression *, expression *);


// return an index to the variable's argument that is better fixed
// in a branching rule for solving a nonconvexity gap
expression *exprMul::getFixVar () 
{return getFixVarBinFun (arglist_ [0], arglist_ [1]);}


/// implied bound processing for expression w = x*y, upon change in
/// lower- and/or upper bound of w, whose index is wind

bool exprMul::impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg) {

  bool res = false;
  int ind;

  if ((arglist_ [ind=0] -> Type () == CONST) || 
      (arglist_ [ind=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [ind] -> Value ();

    ind = arglist_ [1-ind] -> Index ();

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

    return res;
  } else {


    return false; ///////////////////////////////////////////////////////////////

    int xi = arglist_ [0] -> Index (),
        yi = arglist_ [1] -> Index ();

    CouNumber *xl = l + xi,
              *xu = u + xi,
              *yl = l + yi,
              *yu = u + yi,
               wl = l [wind],
               wu = u [wind];

    // w's lower bound 

    if (wl < 0) {

    } else if (wl > 0) {

    }

    // w's upper bound 

  }
}
