/*
 * Name:    CouenneThreeWayBranchObj.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneBranchingObject.hpp>
#include <CouenneThreeWayBranchObj.hpp>


/** \brief Constructor. Get a variable as an argument and set value_
           through a call to operator () of that exprAux.
*/

CouenneThreeWayBranchObj::CouenneThreeWayBranchObj (expression *var, 
						    CouNumber x, 
						    CouNumber l, 
						    CouNumber u): 
  reference_      (var) {

  value_          = x;
  numberBranches_ = 3;

  // Depending on where x, l, and u are, divide bound interval into
  // three and set lcrop_ and rcrop_ accordingly.

  // if l and u are unbounded, crop around x using COUENNE_CROP

  if ((x-l > COUENNE_LARGE_INTERVAL) && 
      (u-x > COUENNE_LARGE_INTERVAL)) {
    lcrop_ = x - COUENNE_CROP;
    rcrop_ = x + COUENNE_CROP;
    first_ = 1;
  }
  else
    if ((x-l > COUENNE_LARGE_INTERVAL) && 
	(u-x < COUENNE_NEAR_BOUND)) {
      lcrop_ = x - COUENNE_CROP;
      rcrop_ = x - COUENNE_LCROP;
      first_ = 2;
    }
    else
      if ((x-l < COUENNE_NEAR_BOUND) && 
	  (u-x > COUENNE_LARGE_INTERVAL)) {
	lcrop_ = x + COUENNE_CROP;
	rcrop_ = x + COUENNE_LCROP;
	first_ = 0;
      }
}


/** \brief Execute the actions required to branch, as specified by the
	   current state of the branching object, and advance the
	   object's state.
	   Returns change in guessed objective on next branch
*/

double CouenneThreeWayBranchObj::branch (OsiSolverInterface * solver) {

  //       -1 if "<= a"  node
  // way =  0 if "[a,b]" node
  //        1 if ">= b"  node

  int way, ind = reference_ -> Index ();

  switch (branchIndex_) {
  case 0: 
    way = first_;   // if first offspring, let first_ decide who's first
    break;
  case 1:
    way = (first_ == 0) ? 1 : 0;
    break;
  case 2:
    way = (first_ == 2) ? 1 : 2;
    break;
  default: 
    printf ("Warning, branchIndex_ has a strange value (%d)\n", branchIndex_);
  }

  way --; // from {0,1,2} to {-1,0,1}

  // set lower or upper bound (round if this variable is integer)

  bool intvar = reference_ -> isInteger ();

  switch (way) {

  case -1: // left interval
    //    printf ("Left branch: x%d <= %.3f\n", ind, lcrop_);
    solver -> setColUpper (ind, intvar ? floor (lcrop_) : lcrop_);
    break;
  case  0: // central interval
    //    printf ("Central branch: %.3f <= x%d <= %.3f\n", lcrop_, ind, rcrop_);
    solver -> setColLower (ind, intvar ? ceil  (lcrop_) : lcrop_);
    solver -> setColUpper (ind, intvar ? floor (rcrop_) : rcrop_);
    break;
  case  1: // right interval
    //    printf ("Right branch: x%d >= %.3f\n", ind, rcrop_);
    solver -> setColLower (ind, intvar ? ceil  (rcrop_) : rcrop_);
    break;
  default:
    printf ("Couenne: branching on nonsense way %d\n", way);
  }

  // TODO: apply bound tightening 

  /*
  if (0) {

  printf (" --> [%.6e,%.6e]\n", l, u);
  printf ("### Branch: x%d %c= %.12f\n", 
  reference_ -> Index (), way ? '>' : '<', value_);

  for (int i=0; i < solver -> getNumCols(); i++)
  printf (" %3d [%.3e,%.3e]\n", i, solver -> getColLower () [i],
  solver -> getColUpper () [i]);
  }*/

  branchIndex_++;

  return 0.; // estimated change in objective function
}
