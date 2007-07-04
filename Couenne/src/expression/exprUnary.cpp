/*
 * Name:    exprUnary.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the unary expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <string>
#include <string.h>

#include <CouenneCutGenerator.hpp>
#include <CouenneProblem.hpp>

#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprAux.hpp>
#include <exprOp.hpp>
#include <exprUnary.hpp>
#include <exprVar.hpp>
#include <exprIVar.hpp>
#include <exprBound.hpp>


// print unary expression


void exprUnary::print (std::ostream &out, 
		       bool descend, 
		       CouenneProblem *p) const {

  if (printPos () == PRE)  out << printOp ();
  out << "("; 
  argument_ -> print (out, descend, p); 
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

exprAux *exprUnary::standardize (CouenneProblem *p) {

  exprAux *subst;

  if ((subst = argument_ -> standardize (p)))
    argument_ = new exprClone (subst);

  return p -> addAuxiliary (this);
}
