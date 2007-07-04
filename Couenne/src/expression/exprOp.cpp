/*
 * Name:    expression.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
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
#include <exprGroup.hpp>

// General N-ary function destructor

exprOp::~exprOp () {

  register expression *elem;

  if (arglist_) {
    for (register int i = nargs_; i--;)
      if ((elem = arglist_ [i]))
	delete elem;

    delete [] arglist_;
    arglist_ = NULL;
  }
}


// print expression

void exprOp::print (std::ostream &out, 
		    bool descend, 
		    CouenneProblem *p) const {
  
  if (printPos () == PRE)
    out << printOp ();

  out << "("; fflush (stdout);
  for (int i=0; i<nargs_; i++) {
    if (arglist_ [i])
      arglist_ [i] -> print (out, descend, p); 
    fflush (stdout);
    if (i < nargs_ - 1) {
      if (printPos () == INSIDE) out << printOp ();
      else                       out << ",";
    }
    fflush (stdout);
  }
  out << ")";
  fflush (stdout);
}


// does this expression depend on variables in varlist?

bool exprOp::dependsOn (int *varlist, int n) {

  if (!varlist) 
    return true;

  for (register int i = nargs_; i--;)
    if (arglist_ [i] -> dependsOn (varlist, n))
      return true;

  return false;
}


/// compare general n-ary expressions

int exprOp::compare (exprOp  &e1) {

  int c0 = code (),
      c1 = e1. code ();

  if (c0 < c1) return -1;
  if (c0 > c1) return  1;

  // have to compare arguments one by one
  if (nargs_ < e1.nargs_) return -1;
  if (nargs_ > e1.nargs_) return  1;

  // not an exprGroup, compare arguments
  for (register int i = nargs_; i--;) {

    int res = arglist_ [i] -> compare (*(e1. ArgList () [i]));
    if (res) return res;
  }

  // last chance, this might be an exprGroup
  if (c0==COU_EXPRGROUP) {

    exprGroup *ne0 = dynamic_cast <exprGroup *> (this),
              *ne1 = dynamic_cast <exprGroup *> (&e1);

    return ne0 -> compare (*ne1);
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprOp::rank (CouenneProblem *p) {

  int maxrank = -1;

  for (register expression **al = arglist_ + nargs_; 
       al-- > arglist_;) {
    register int r = (*al) -> rank (p);
    if (r > maxrank) maxrank = r;
  }

  return (1 + maxrank);
}


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)
//
// For the base exprOp class we only do the first part (for argument
// list components only), and the calling class (Sum, Sub, Mul, Pow,
// and the like) will do the part for its own object

exprAux *exprOp::standardize (CouenneProblem *p) {

  exprVar *subst;

  for (register int i = nargs_; i--;)

    if ((subst = arglist_ [i] -> standardize (p)))
      arglist_ [i] = new exprClone (subst);

  return NULL;
}
