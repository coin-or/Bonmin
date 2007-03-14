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

CouNumber CouenneBranchingObject::alpha_ = 0.1;


/** \brief Constructor. Get a variable as an argument and set value_
           through a call to operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (expression *var): 
  reference_ (var) {

  // set the branching value at the current point 
  //  value_ = expression::Variable (reference_ -> Index ());

  int index = reference_ -> Index ();

  CouNumber 
    x = expression::Variable (index),   // current solution
    l = expression::Lbound   (index),   //         lower bound
    u = expression::Ubound   (index),   //         upper
    alpha = CouenneBranchingObject::Alpha ();

  if ((x > l + COUENNE_EPS) && 
      (x < u - COUENNE_EPS))
    // infinite (at least on one side) bound interval, but x is not
    // at the boundary
    value_ = x;
  else // current point is at one of the bounds
    if ((l > - COUENNE_INFINITY + 1) &&
	(u <   COUENNE_INFINITY - 1)) 
      // finite bounds, apply midpoint rule
      value_ = alpha * x + (1-alpha) * (l+u) / 2.;
    else 
      // infinite bound interval, x is at the boundary
      // push it inwards
      // TODO: look for a proper value for the displacement 
      if (fabs (x-l) < COUENNE_EPS) value_ += (1+fabs (l)) / 2.; 
      else                          value_ -= (1+fabs (u)) / 2.; 

  if (0) {
    printf ("Branch::constructor: ");
    reference_ -> print (std::cout);
    printf (" on %f [%f,%f]\n", 
	    value_, 
	    expression::Lbound (reference_ -> Index ()),
	    expression::Ubound (reference_ -> Index ()));
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

  int way = (!branchIndex_) ? firstBranch_ : !firstBranch_;
  /*
  if (way) {

    if (value_ < solver -> getColLower () [reference_ -> Index ()])
      printf ("Nonsense   up-branch: [ %f ,(%f)] -> %f\n", 
	      solver -> getColLower () [reference_ -> Index()], 
	      solver -> getColUpper () [reference_ -> Index()], value_);
    else
      if (value_ < solver -> getColLower () [reference_ -> Index ()] + COUENNE_EPS)
	printf ("## WEAK    up-branch: [ %f ,(%f)] -> %f\n", 
		solver -> getColLower () [reference_ -> Index()], 
		solver -> getColUpper () [reference_ -> Index()], value_);
  } else {

    if (value_ > solver -> getColUpper () [reference_ -> Index()])
      printf ("Nonsense down-branch: [(%f), %f ] -> %f\n", 
	      solver -> getColLower () [reference_ -> Index()], 
	      solver -> getColUpper () [reference_ -> Index()], value_);
    else
      if (value_ > solver -> getColUpper () [reference_ -> Index()] + COUENNE_EPS)
	printf ("## WEAK  down-branch: [(%f), %f ] -> %f\n", 
		solver -> getColLower () [reference_ -> Index()], 
		solver -> getColUpper () [reference_ -> Index()], value_);
  }
  */
  if (way) // ">=" node, set lower bound (round if this variable is integer)
    solver -> setColLower (reference_ -> Index (), 
			   reference_ -> isInteger () ? ceil  (value_) : value_);
  else     // "<=" node, set upper bound (ditto)
    solver -> setColUpper (reference_ -> Index (), 
			   reference_ -> isInteger () ? floor (value_) : value_);

  //  printf ("################################# Branch: x%d %c= %.12f\n", 
  //	  reference_ -> Index (), way ? '>' : '<', value_);

  branchIndex_++;
  return 0.; // estimated change in objective function
}
