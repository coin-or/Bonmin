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

void exprOp::print (std::ostream      &out = std::cout, 
		    const std::string &op  = "unknown", 
		    enum pos           pos = PRE)        const 
{
  if (pos == PRE)
    out << op;

  out << "("; fflush (stdout);
  for (int i=0; i<nargs_; i++) {
    if (arglist_ [i])
      arglist_ [i] -> print (out); 
    fflush (stdout);
    if (i < nargs_ - 1) {
      if (pos == INSIDE) out << op;
      else               out << ",";
    }
    fflush (stdout);
  }
  out << ")";
  fflush (stdout);
}


// does this expression depend on variables in varlist?

bool exprOp::dependsOn (int *varlist = NULL, int n = 1) {

  for (register int i = nargs_; i--;)
    if (arglist_ [i] -> dependsOn (varlist, n))
      return true;

  return false;
}

///
int exprOp::compare (exprOp  &e1) {

  int c0 = code (),
      c1 = e1. code ();

  if      (c0 < c1) return -1;
  else if (c0 > c1) return  1;
  else { // have to compare arguments one by one

    if      (nargs_ < e1.nargs_) return -1;
    else if (nargs_ > e1.nargs_) return 1;

    for (register int i = nargs_; i--;) {

      int res = arglist_ [i] -> compare (*(e1. ArgList () [i]));
      if (res) return res;
    }
    return 0;
  }
}
