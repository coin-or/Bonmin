/*
 * Name:    branchExprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for qudratic forms
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>


//#define LARGE_BOUND 9.9e12
//#define BOUND_WING 1

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprQuad::selectBranch (expression *w, 
				  const OsiBranchingInformation *info,
				  int &ind, 
				  double * &brpts, 
				  int &way) {

  return 0;

  const CouNumber *x = info -> solution_, 
                  *l = info -> lower_,
                  *u = info -> upper_; 




#if 0

  int xi = arglist_ [0] -> Index (),
      yi = arglist_ [1] -> Index ();

  if ((xi < 0) || (yi < 0)) {
    printf ("Couenne, exprMul::selectBranch: arguments of exprMul have negative index\n");
    exit (-1);
  }

  CouNumber x0 = expression::Variable (xi),   y0 = expression::Variable (yi), 
            xl = expression::Lbound   (xi),   yl = expression::Lbound   (yi),  
            xu = expression::Ubound   (xi),   yu = expression::Ubound   (yi);  

  // First, try to avoid infinite bounds for multiplications, which
  // make them pretty hard to deal with

  if ((xl < - COUENNE_INFINITY) && 
      (xu >   COUENNE_INFINITY)) { // first

    ind = xi;

    // branch around current point. If it is also at a crazy value,
    // reset it close to zero.

    brpts = (double *) malloc (2 * sizeof (double));
    CouNumber curr = x0;

    if (fabs (curr) >= LARGE_BOUND) curr = 0;

    brpts [0] = - fabs (curr) - BOUND_WING;
    brpts [1] =   fabs (curr) + BOUND_WING;

    way = THREE_CENTER;

    return - COUENNE_INFINITY; // tell caller not to set infeasibility to this
  }

  if ((yl < - COUENNE_INFINITY) && 
      (yu >   COUENNE_INFINITY)) { // and second factor

    ind = yi;

    // branch around current point. If it is also at a crazy value,
    // reset it close to zero.

    brpts = (double *) malloc (2 * sizeof (double));
    CouNumber curr = y0;

    if (fabs (curr) >= LARGE_BOUND) curr = 0;

    brpts [0] = - fabs (curr) - BOUND_WING;
    brpts [1] =   fabs (curr) + BOUND_WING;

    way = THREE_CENTER;

    return - COUENNE_INFINITY; // tell caller not to set infeasibility to this
  }

  // there is at least one finite corner of the bounding box

  ind = -1;
  return 0.;
#endif
}
