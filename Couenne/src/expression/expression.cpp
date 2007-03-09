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

  OsiRowCut *cut;

  if ((cut = cg -> createCut (currValue_, 0, w -> Index (), 1.)))
    cs.insert (cut);
}


/// compare generic expression with other expression
int expression::compare (expression &e1) {

  if      (code () >= COU_EXPRUNARY) 
    if (e1.code () >= COU_EXPRUNARY) {

      exprUnary *ne0 = dynamic_cast <exprUnary *> (this);
      exprUnary *ne1 = dynamic_cast <exprUnary *> (&e1);
      return ne0 -> compare (*ne1);
    }
    else return 1;
  else if   (code () >= COU_EXPROP)
    if   (e1.code () >= COU_EXPROP)
      if (e1.code () >= COU_EXPRUNARY) return -1;
      else {

      exprOp *ne0 = dynamic_cast <exprOp *> (this);
      exprOp *ne1 = dynamic_cast <exprOp *> (&e1);
      return ne0 -> compare (*ne1);
      }
    else if (e1.code () >= COU_EXPROP) return -1;
    else ;
  else ;

  int c0 = code (),
      c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;
  else { // it is either a constant or a variable

    int i1 = e1. Index ();

    if      (Index () < i1) return -1;
    else if (Index () > i1) return  1;
    else if (i1 != -1)      return  0; // same variables
    else { // both are constants

      CouNumber v1 = e1. Value ();

      if      (currValue_ < v1) return -1;
      else if (currValue_ > v1) return  1;
      else                      return  0;
    }
  }
}


///
int expression::compare (exprCopy &c)
{return compare (const_cast <expression &> (*(c. Original ())));}
