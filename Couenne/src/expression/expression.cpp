/*
 * Name:    expression.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneCutGenerator.hpp>
#include <CouenneProblem.hpp>

#include <CouenneTypes.hpp>
#include <expression.hpp>
#include <exprAux.hpp>
#include <exprOp.hpp>
#include <exprUnary.hpp>


// static vectors for evaluation, see their description in
// expression.hpp

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
			      OsiCuts &cs, const CouenneCutGenerator *cg, 
			      t_chg_bounds *chg, int,
			      CouNumber, CouNumber) {

  if (cg -> isFirst ())
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


/// replace occurrence of a variable with another variable
void exprCopy::replace (exprVar *orig, exprVar *aux) {

  if ((copy_ -> Index () == orig -> Index ()) &&
      (copy_ -> Type  () != AUX)) {

    delete copy_;
    copy_ = new exprClone (aux);
  }
}


/// dependence on variable set: return cardinality of subset of the
/// set of indices in first argument which occur in expression. 
int expression::dependsOn (int *ind, int n, 
			   CouenneProblem *p, 
			   enum dig_type type) {

  std::set <int> 
    indlist (ind, ind + n), 
    deplist,
    intersectn;

  DepList (deplist, type, p);

  std::set_intersection (indlist .begin (), indlist .end (), 
			 deplist .begin (), deplist .end (),
			 std::inserter (intersectn, intersectn.begin ()));

  return intersectn.size();
}
