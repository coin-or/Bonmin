/*
 * Name:    simplify.cpp
 * Author:  Pietro Belotti
 * Purpose: symbolic expression simplifier
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "expression.hpp"
#include "exprOp.hpp"
#include "exprUnary.hpp"
#include "exprConst.hpp"


// simplify n-ary expression f (g_1(x), g_2(x)... g_n(x))

expression *exprOp:: simplify () {

  //  Simplify arguments g_1(x), g_2(x)... g_n(x) first
  for (int i=0; i<nargs_; i++) {

    expression *subst;

    if ((subst = (arglist_ [i]) -> simplify ())) {

      delete arglist_ [i];
      arglist_ [i] = subst;
    }
  }

  return NULL;
}


// simplify unary operators

expression *exprUnary:: simplify () {

  register expression *subst;

  // Simplify argument g(x) of this expression f(g(x))
  if ((subst = argument_ -> simplify ())) {

    delete argument_;
    argument_ = subst;

    // g(x) is a constant k, therefore return f (k)
    if (subst -> Type () == CONST) {

      expression *ret = new exprConst (operator () ());
      argument_ = NULL;
      delete subst;

      return ret;
    } else return NULL;
  } else return NULL;
}
