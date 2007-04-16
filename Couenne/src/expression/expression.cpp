/*
 * Name:    expression.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <string>
#include <string.h>

#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprIVar.h>
#include <exprBound.h>

// static vectors for evaluation, see their description in
// expression.h

CouNumber  expression::stack [STACK_SIZE];
CouNumber *expression::sp = stack;

CouNumber *expression::variables_ = NULL;
CouNumber *expression::lbounds_   = NULL;
CouNumber *expression::ubounds_   = NULL;


// Get lower and upper bound of a generic expression

void expression::getBounds (expression *&lb, expression *&ub) {

  lb = new exprConst (- COUENNE_INFINITY);
  ub = new exprConst (  COUENNE_INFINITY);
}


// generate one cut for a constant

void exprConst::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg) {

  cg -> createCut (cs, currValue_, 0, w -> Index (), 1.);
}


/// compare generic expression with other expression
int expression::compare (expression &e1) {

  register int c0 = code (),
               c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;

  // same code, check arguments

  if (c0 >= COU_EXPRUNARY) { // both are exprUnary's

    exprUnary *ne0 = dynamic_cast <exprUnary *> (this);
    exprUnary *ne1 = dynamic_cast <exprUnary *> (&e1);

    return ne0 -> compare (*ne1);
  }

  if (c0 >= COU_EXPROP) { // both are exprOp's

    exprOp *ne0 = dynamic_cast <exprOp *> (this);
    exprOp *ne1 = dynamic_cast <exprOp *> (&e1);

    return ne0 -> compare (*ne1);
  }

  // expressions are both variables or constants

  register int i0 =     Index (),
               i1 = e1. Index ();

  if (i0 < i1) return -1;
  if (i0 > i1) return  1;
  if (i0 >= 0) return  0; // same variables

  // both are numbers
  register CouNumber v1 = e1. Value ();

  if (currValue_ < v1) return -1;
  if (currValue_ > v1) return  1;

  return  0;
}


/// compare expressions (used in bsearch within CouenneProblem::standardize)
int expression::compare (exprCopy &c)
{return compare (const_cast <expression &> (*(c. Original ())));}
