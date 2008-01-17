/*
 * Name:    exprUnary.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the unary expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneTypes.hpp"
#include "exprUnary.hpp"
#include "exprVar.hpp"


// print unary expression
void exprUnary::print (std::ostream &out, 
		       bool descend) const {

  if (printPos () == PRE)  out << printOp ();
  out << "("; 
  argument_ -> print (out, descend); 
  out << ")";
  if (printPos () == POST) out << printOp ();
}


/// comparison when looking for duplicates
int exprUnary::compare (exprUnary  &e1) { 

  int c0 = code (),
      c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;
  else // have to compare arguments 
    return argument_ -> compare (*(e1.argument_));
}


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)
exprAux *exprUnary::standardize (CouenneProblem *p, bool addAux) {

  exprAux *subst;

  if ((subst = argument_ -> standardize (p)))
    argument_ = new exprClone (subst);

  return (addAux ? (p -> addAuxiliary (this)) : new exprAux (this));
}


/// replace variable with other
void exprUnary::replace (exprVar *x, exprVar *w) {

  if (argument_ -> Type () == VAR) {
    if (argument_ -> Index () == x -> Index ()) {
      delete argument_;
      argument_ = new exprClone (w);
    }
  } else argument_ -> replace (x, w);
}


/// is this expression integer?
bool exprUnary::isInteger () {

  // only check if argument is, *at this point in the algorithm*,
  // constant -- due to branching rules, for instance. If so, check if
  // the corresponding evaluated expression is integer.

  CouNumber al, au;
  argument_ -> getBounds (al, au);

  if (fabs (al - au) < COUENNE_EPS) { // argument is constant

    register CouNumber fval = (F ()) (al); 

    // check if f(lb=ub) is integer
    if (fabs (COUENNE_round (fval) - fval) < COUENNE_EPS)
      return true;
  }

  return false;
}
