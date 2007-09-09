/*
 * Name:    exprOp.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for multivariate function class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <string>
#include <string.h>

#include <CouenneProblem.hpp>

#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprAux.hpp>
#include <exprOp.hpp>
#include <exprUnary.hpp>
#include <exprGroup.hpp>
#include <exprQuad.hpp>

// General N-ary function destructor

exprOp::~exprOp () {

  if (arglist_) {
    for (register int i = nargs_; i--;)
      if (arglist_ [i])
	delete arglist_ [i];

    delete [] arglist_;
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

int exprOp::dependsOn (int *varlist, int n) {

  int tot = 0;

  for (register int i = nargs_; i--;)
    tot += arglist_ [i] -> dependsOn (varlist, n);

  return tot;
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

  // last chance, this might be an exprGroup or derived
  if ((c0 == COU_EXPRGROUP) ||
      (c0 == COU_EXPRQUAD)) {

    exprGroup *ne0 = dynamic_cast <exprGroup *> (this),
              *ne1 = dynamic_cast <exprGroup *> (&e1);

    int cg = ne0 -> compare (*ne1);

    if (cg) return cg; // exprGroup

    // last chance, the two are quadratic forms

    if (c0 == COU_EXPRQUAD) {

      exprQuad *ne0 = dynamic_cast <exprQuad *> (this),
   	       *ne1 = dynamic_cast <exprQuad *> (&e1);

      return ne0 -> compare (*ne1);
    }
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

  return (maxrank);
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

exprAux *exprOp::standardize (register CouenneProblem *p) {

  register exprVar *subst;

  for (register int i = nargs_; i--;)
    if ((subst = arglist_ [i] -> standardize (p)))
      arglist_ [i] = new exprClone (subst);

  return NULL;
}


/// replace variable x with new (aux) w
void exprOp::replace (exprVar *x, exprVar *w) {

  /*printf ("replacing "); fflush (stdout);
  x -> print (); printf (" with "); fflush (stdout);
  w -> print (); printf (" in "); fflush (stdout);
  print (); fflush (stdout); printf ("\n");*/

  register expression **al = arglist_;

  for (register int i=nargs_; i--; al++)
    if ((*al) -> Type () == VAR) {
      if ((*al) -> Index () == x -> Index ()) {
	delete *al;
	*al = new exprClone (w);
      }
    } else (*al) -> replace (x, w);
}
