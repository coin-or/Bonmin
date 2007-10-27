/*
 * Name:    exprOp.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for multivariate function class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.hpp>

#include <expression.hpp>
#include <exprAux.hpp>
#include <exprOp.hpp>
#include <exprGroup.hpp>
#include <exprQuad.hpp>


// General N-ary function destructor

exprOp::~exprOp () {

  if (arglist_) {

    for (int i = nargs_; i--;) {

      expression *&arg = arglist_ [i];

      if (arg) {
	delete arg;
	arg = NULL;
      }
    }

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

exprAux *exprOp::standardize (register CouenneProblem *p, bool addAux) {

  register exprVar *subst;

  for (register int i = nargs_; i--;)
    if ((subst = arglist_ [i] -> standardize (p)))
      arglist_ [i] = new exprClone (subst);

  return NULL;
}


/// replace variable x with new (aux) w
void exprOp::replace (exprVar *x, exprVar *w) {

  expression **al = arglist_;
  int index = x -> Index ();

  for (register int i = nargs_; i--; al++)

    switch ((*al) -> Type () == VAR) {

    case VAR:
      if ((*al) -> Index () == index) {
	delete *al;
	*al = new exprClone (w);
      }
      break;

    case UNARY:
    case N_ARY:
      (*al) -> replace (x, w);
      break;

    default:
      break;
    }
}


/// is this expression integer?
bool exprOp::isInteger () {

  for (int i = nargs_; i--;)

    if (!(arglist_ [i] -> isInteger ())) { // this argument is not integer

      // last chance: check  if constant and integer

      expression *lb, *ub;

      arglist_ [i] -> getBounds (lb, ub);
      CouNumber lv = (*lb) ();

      if ((fabs (lv - (*ub) ()) > COUENNE_EPS) ||
	  (fabs (COUENNE_round (lv) - lv) < COUENNE_EPS))
      return false;
    }

  return true;
}
