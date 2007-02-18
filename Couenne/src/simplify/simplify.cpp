/*
 * Name:    simplify.cpp
 * Author:  Pietro Belotti
 * Purpose: symbolic expression simplifier
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <expression.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprConst.h>


// simplify n-ary expression f (g_1(x), g_2(x)... g_n(x))

expression *exprOp:: simplify () {

  //  print (std::cout, "op", PRE); printf (" --> ");

  //  Simplify arguments g_1(x), g_2(x)... g_n(x) first
  for (register int i=0; i<nargs_; i++) {

    register expression *subst;

    if ((subst = (arglist_ [i]) -> simplify ())) {

      delete arglist_ [i];
      arglist_ [i] = subst;
    }
  }

  //  print (std::cout, "op", PRE); printf ("\n");

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
    }
    else return NULL;
  }
  else return NULL;
}
