/*
 * Name:    CouenneBranchingObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CoinHelperFunctions.hpp>

#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>

//#define DEBUG

/// make branching point $\alpha$ away from current point:
/// bp = alpha * current + (1-alpha) * midpoint

CouNumber CouenneBranchingObject::alpha_ = 0.25;


/** \brief Constructor. 
 *
 * Get a variable as an argument and set value_ through a call to
 * operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (int index, int way, CouNumber brpoint, bool isint): 

  index_   (index),
  integer_ (isint) {

  firstBranch_ =  (way == TWO_LEFT)      ? 0 : 
                 ((way == TWO_RIGHT)     ? 1 : 
                 ((CoinDrand48 () < 0.5) ? 0 : 1));

  if (index_ < 0) {
    printf ("Couenne: CouenneBranchingObject has negative reference's index\n");
    exit (-1);
  }

  CouNumber
    x = expression::Variable (index_),   // current solution
    l = expression::Lbound   (index_),   //         lower bound
    u = expression::Ubound   (index_),   //         upper
    alpha = CouenneBranchingObject::Alpha ();

  if (fabs (brpoint) < COUENNE_INFINITY) 
    x = brpoint;

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

#ifdef DEBUG
  if (fabs (u-l) < COUENNE_EPS)
    printf ("#### Warning, interval is really tiny\n");
#endif

  if ((x-l < COUENNE_NEAR_BOUND) ||
      (u-x < COUENNE_NEAR_BOUND))
    if (u < COUENNE_INFINITY)
      if (l > - COUENNE_INFINITY)    value_ = alpha * x + (1. - alpha) * (l + u) / 2.;
      else                           value_ = (u<0) ? (u-1) : u/2;
    else if (l > - COUENNE_INFINITY) value_ = (l>0) ? (l+1) : l/2;
    else                             value_ = 0;
  else value_ = x;

#ifdef DEBUG
  printf ("=== x%d will branch on %g (at %g) [%g,%g]\n", 
	  index_, value_, 
	  expression::Variable (index_),
	  expression::Lbound   (index_),
	  expression::Ubound   (index_));
#endif
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

#ifdef DEBUG

  CouNumber l = solver -> getColLower () [index_],
            u = solver -> getColUpper () [index_];

  if (way) {
    if      (value_ < l)             printf ("Nonsense up-br: [ %.8f ,(%.8f)] -> %.8f\n", l,u,value_);
    else if (value_ < l+COUENNE_EPS) printf ("## WEAK  up-br: [ %.8f ,(%.8f)] -> %.8f\n", l,u,value_);
  } else {
    if      (value_ > u)             printf ("Nonsense dn-br: [(%.8f), %.8f ] -> %.8f\n", l,u,value_);
    else if (value_ > u+COUENNE_EPS) printf ("## WEAK  dn-br: [(%.8f), %.8f ] -> %.8f\n", l,u,value_);
  }
#endif

  if (!way) solver -> setColUpper (index_, integer_ ? floor (value_) : value_); // down branch
  else      solver -> setColLower (index_, integer_ ? ceil  (value_) : value_); // up   branch

  // TODO: apply bound tightening to evaluate change in dual bound

  //    printf (" --> [%.6e,%.6e]\n", l, u);

#ifdef DEBUG
  printf ("### Branch: x%d %c= %g\n", 
  	  index_, way ? '>' : '<', value_);
#endif

  branchIndex_++;
  return 0.; // estimated change in objective function
}
