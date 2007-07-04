/*
 * Name:    simplifiers.cpp
 * Author:  Pietro Belotti
 * Purpose: simplifiers for main operators (+,*,-)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprOp.hpp>
#include <exprConst.hpp>
#include <CouennePrecisions.h>


//
// shrink argument list
//
// used by + and * (for now), accepts a constant resulting from
// applying an operator to the constants in the list of (pointers to)
// function arguments contained in el. The constant is inserted in the
// list if the result is not equal to null_element or if there are
// other non-constant terms in the arglist.
// 
// Example: f(x) + 3 + g(x) + 2 + 4 
//
// el    = {pf, NULL, pg, NULL, NULL}
// nargs = 5 
// c     = 3 + 2 + 4 = 9
// null_element = 0 (for sums)
// 
// where pf and pg are pointers to expression containing f and g,
// resp.
//
// Result: el and nargs are changed to
//
// el    = {pf, p9, pg}
// nargs = 3 
//
// Another example: f(x) + 2 + g(x) + (-4) + 2
// Result:
// el    = {pf, pg}
// nargs = 2
//
// Another example: f(x) * 3 * g(x) * 2 
//
// el    = {pf, NULL, pg, NULL}
// nargs = 4
// c     = 3 * 2 = 6 != null_element = 1 (for multiplications)
// Result:
// el    = {pf, p6, pg}
// nargs = 3
//

int exprOp::shrink_arglist (CouNumber c, CouNumber null_element) {

  register int i=0, j=0;

  bool one_fun = false;

  // find first NULL spot (left by some constant)
  while ((i < nargs_) && (arglist_ [i])) 
    i++; 

  // no spots, leave
  if (i==nargs_) 
    return 0;

  // check if there is at least one non-constant expression
  for (register int k=nargs_; k--;) 
    if (arglist_ [k]) {
      one_fun = true;
      break;
    }

  // add constant term if c is not null w.r.t the operation or if it
  // would be an empty operand list otherwise
  if ((fabs (c - null_element) > COUENNE_EPS) || !one_fun)
    arglist_ [i++] = new exprConst (c);

  j = i;

  //  printf ("////////////////////////////// j=%d\n", j);

  // now shift back all operands to compress argument list
  while (i < nargs_) {

    while ((i < nargs_) && !(arglist_ [i])) 
      i++;

    if (i < nargs_) 
      one_fun = true;

    while ((i < nargs_) && (arglist_ [i]))
      arglist_ [j++] = arglist_ [i++]; 
  }

  nargs_ = j;

  // only say shrinking simplified arg list if there is just the
  // constant
  return (nargs_ == 1);// && ((fabs (c - null_element) > COUENNE_EPS) || !one_fun));
}
