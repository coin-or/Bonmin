/*
 * Name:    getFixVarBinFun.cpp
 * Author:  Pietro Belotti
 * Purpose: return which argument, in a binary function
 *          (multiplications and divisions), should be 
 *          branched on.
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "expression.hpp"
#include "domain.hpp"

// return an index to the variable's argument that is better fixed
// in a branching rule for solving a nonconvexity gap

expression *getFixVarBinFun (expression *arg0, expression *arg1) {

  if (arg0 -> Type () == CONST) return arg1; // Case c*y or c/y
  if (arg1 -> Type () == CONST) return arg0; // Case x*d or x/d

  // Case x*y, x/y: return variable closest to midpoint in bound interval

  // get variables' indices

  int 
    index0 = arg0 -> Index (),
    index1 = arg1 -> Index ();

  DomainPoint *point = arg0 -> domain () -> current ();

  // get bounds of both variables

  CouNumber 
    l0 = point -> lb (index0), l1 = point -> lb (index1), 
    u0 = point -> ub (index0), u1 = point -> ub (index1), 
    x0 = point -> x  (index0), x1 = point -> x (index1), 
    //l0 = expression::Lbound   (index0), l1 = expression::Lbound   (index1), 
    //u0 = expression::Ubound   (index0), u1 = expression::Ubound   (index1),
    //x0 = expression::Variable (index0), x1 = expression::Variable (index1),

    delta0 = u0-l0, 
    delta1 = u1-l1;

  if        (fabs (delta0) < COUENNE_EPS) {
    if      (fabs (delta1) > COUENNE_EPS) return arg1;
  } else if (fabs (delta1) < COUENNE_EPS) return arg0;

  // We have a full-dimensional (possibly unlimited) bounding box,
  // or a very tiny square

  { // branch on variable far from bounds if other is near bounds

    register bool 
      xl0 = (fabs (x0-l0) < COUENNE_EPS), 
      xu0 = (fabs (x0-u0) < COUENNE_EPS),
      xl1 = (fabs (x1-l1) < COUENNE_EPS), 
      xu1 = (fabs (x1-u1) < COUENNE_EPS);

    if      ((xl0 || xu0) && !(xl1 || xu1)) return arg1;
    else if ((xl1 || xu1) && !(xl0 || xu0)) return arg0;
  }

  // The bounding box may be unlimited but at least no variable is
  // close to the bound

  if ((  l0 > - COUENNE_INFINITY) && (u0 < COUENNE_INFINITY)) {
    if ((l1 > - COUENNE_INFINITY) && (u1 < COUENNE_INFINITY)) {

      // Ideal situation, bounding box is finite.  Branch on the one
      // closest to the midpoint

      CouNumber 
	normdist0 = fabs (0.5 * (u0+l0) - x0) / delta0,
	normdist1 = fabs (0.5 * (u1+l1) - x1) / delta1;

      if (normdist0 < normdist1) 
	return arg0;
      else return arg1;

    } else return arg0; // y has at least one infinite bound
  }
  else {

    if ((l1 > - COUENNE_INFINITY) && 
	(u1 <   COUENNE_INFINITY))
      return arg1;
    else {

      // Both bound intervals are unlimited, we choose the variable
      // with the "most internal" point

      CouNumber distance0,
	distance1 = distance0 = - COUENNE_INFINITY;

      if      (l0 > - COUENNE_INFINITY) distance0 = x0 - l0;
      else if (u0 <   COUENNE_INFINITY) distance0 = u0 - x0;
	
      if      (l1 > - COUENNE_INFINITY) distance1 = x1 - l1;
      else if (u1 <   COUENNE_INFINITY) distance1 = u1 - x1;

      if (distance0 > distance1) 
	return arg0;
      else return arg1;
    }
  }

}
