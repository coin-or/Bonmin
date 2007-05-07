/*
 * Name:    CouenneBranchingObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneBranchingObject.hpp>

/// make branching point $\alpha$ away from current point:
/// bp = alpha * current + (1-alpha) * midpoint

CouNumber CouenneBranchingObject::alpha_ = 0.25;


/** \brief Constructor. Get a variable as an argument and set value_
           through a call to operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (expression *var): 

  reference_ (var) {

  // set the branching value at the current point 
  //  value_ = expression::Variable (reference_ -> Index ());

  int index = reference_ -> Index ();

  if (index < 0)
    printf ("Couenne: Warning, CouenneBranchingObject has negative reference's index\n");

  long double
    x = expression::Variable (index),   // current solution
    l = expression::Lbound   (index),   //         lower bound
    u = expression::Ubound   (index),   //         upper
    alpha = CouenneBranchingObject::Alpha ();

  if      (x<l) x = l;
  else if (x>u) x = u;

  // This two-way branching rule is only applied when both lower and
  // upper bound are finite. Otherwise, a CouenneThreeWayBranchObj is
  // used (see CouenneThreeWayBranchObj.hpp).
  //
  // The rule is as follows:
  //
  // - if x is well inside the interval (both bounds are infinite or
  // there is a difference of at least COU_NEAR_BOUND), set
  // value_ to x;
  //
  // - otherwise, try to get far from bounds by setting value_ to a
  // convex combination of current and midpoint
  //
  // TODO: consider branching value that maximizes distance from
  // current point (how?)

  if (fabs (u-l) < COUENNE_EPS)
    printf ("Warning, interval is really tiny\n");

  if ((x-l < COUENNE_NEAR_BOUND) ||
      (u-x < COUENNE_NEAR_BOUND))
    value_ = alpha * x + (1. - alpha) * (l + u) / 2.;
  else value_ = x;

  /*
  if ((x > l + COUENNE_NEAR_BOUND) && 
      (x < u - COUENNE_NEAR_BOUND))      // x is not at the boundary
    value_ = x;
  else // current point is at one of the bounds
    if ((l > - COUENNE_INFINITY) &&
	(u <   COUENNE_INFINITY))
      // finite bounds, apply midpoint rule
      value_ = alpha * x + (1. - alpha) * (l + u) / 2.;

    else
      // infinite (one direction) bound interval, x is at the boundary
      // push it inwards
      // TODO: look for a proper displacement 
      if (fabs (x-l) < COUENNE_EPS) value_ = l + (1+fabs (l)) / 2.; 
      else                          value_ = u - (1+fabs (u)) / 2.; 
  */

  if (0) {
    printf ("Branch::constructor: ");
    reference_ -> print (std::cout);
    printf (" on %.6e (%.6e) [%.6e,%.6e]\n", 
	    value_, 
	    expression::Variable (reference_ -> Index ()),
	    expression::Lbound   (reference_ -> Index ()),
	    expression::Ubound   (reference_ -> Index ()));
  }
}


/** \brief Execute the actions required to branch, as specified by the
	   current state of the branching object, and advance the
	   object's state.
	   Returns change in guessed objective on next branch
*/

double CouenneBranchingObject::branch (OsiSolverInterface * solver) {

  // way = 0 if "<=" node, 
  //       1 if ">=" node

  int way = (!branchIndex_) ? firstBranch_ : !firstBranch_,
      ind = reference_ -> Index ();

  CouNumber l = solver -> getColLower () [ind],
            u = solver -> getColUpper () [ind];

  if (way) {
    if      (value_ < l)             printf ("Nonsense up-branch: [ %f ,(%f)] -> %f\n", l,u,value_);
    else if (value_ < l+COUENNE_EPS) printf ("## WEAK  up-branch: [ %f ,(%f)] -> %f\n", l,u,value_);
  } else {
    if      (value_ > u)             printf ("Nonsense dn-branch: [(%f), %f ] -> %f\n", l,u,value_);
    else if (value_ > u+COUENNE_EPS) printf ("## WEAK  dn-branch: [(%f), %f ] -> %f\n", l,u,value_);
  }

  if (0) printf (" [%.6e,%.6e]", l, u);

  // ">=" if way=1, "<=" if way=0 set lower or upper bound (round if this variable is integer)
  if (way) solver -> setColLower (ind, reference_ -> isInteger () ? ceil  (value_) : value_);
  else     solver -> setColUpper (ind, reference_ -> isInteger () ? floor (value_) : value_);

  // TODO: apply bound tightening 

  if (0) {

    printf (" --> [%.6e,%.6e]\n", l, u);
    printf ("### Branch: x%d %c= %.12f\n", 
	    reference_ -> Index (), way ? '>' : '<', value_);

    /*for (int i=0; i < solver -> getNumCols(); i++)
      printf (" %3d [%.3e,%.3e]\n", i, solver -> getColLower () [i],
      solver -> getColUpper () [i]);*/
  }

  branchIndex_++;
  return 0.; // estimated change in objective function
}
