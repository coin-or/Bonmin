/*
 * Name:    CouenneBranchingObject.cpp
 * Author:  Pietro Belotti
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneBranchingObject.hpp>


/** \brief Constructor. Get a variable as an argument and set value_
           through a call to operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (expression *var): 
  reference_ (var)
  {value_ = (*reference_) ();} // set the branching value at the current point 


/** \brief Execute the actions required to branch, as specified by the
	   current state of the branching object, and advance the
	   object's state.
	   Returns change in guessed objective on next branch
*/

double CouenneBranchingObject::branch (OsiSolverInterface * solver) {

  // way = 0 if "<=" node, 
  //       1 if ">=" node

  int way = (!branchIndex_) ? firstBranch_ : !firstBranch_;

  if (way) // ">=" node, set lower bound (round if this variable is integer)
    solver -> setColLower (reference_ -> Index(), 
			   reference_ -> isInteger() ? ceil  (value_) : value_);
  else     // "<=" node, set upper bound (ditto)
    solver -> setColUpper (reference_ -> Index(), 
			   reference_ -> isInteger() ? floor (value_) : value_);

  return 0.; // estimated gain in objective function
}
