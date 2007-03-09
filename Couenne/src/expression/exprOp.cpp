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
#include <exprGroup.h>

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

  // last chance, this might be an exprGrouip
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
