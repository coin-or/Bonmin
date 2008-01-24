/*
 * Name:    CouenneBranchingObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

/// make branching point $\alpha$ away from current point:
/// bp = alpha * current + (1-alpha) * midpoint

CouNumber CouenneBranchingObject::alpha_ = 0.25;

#define AGGR_MUL 2


/// computes a not-too-bad point where to branch, in the "middle" of an interval
CouNumber midInterval (CouNumber curr, CouNumber l, CouNumber u) {

  CouNumber x = curr;

  if (u < l + COUENNE_EPS)
    return (0.5 * (l + u));

  if      (x<l) x = l;
  else if (x>u) x = u;

  CouNumber alpha = CouenneBranchingObject::Alpha ();

  if   (l < -COUENNE_INFINITY / 10)
    if (u >  COUENNE_INFINITY / 10) return x; // 0.                                    // ]-inf,+inf[
    else                            return ((x < -COUENNE_EPS) ? (AGGR_MUL * (-1+x)) : // ]-inf,u]
					    (x >  COUENNE_EPS) ? 0. : -AGGR_MUL);
  else
    if (u >  COUENNE_INFINITY / 10) return ((x >  COUENNE_EPS) ? (AGGR_MUL *  (1+x)) : // [l,+inf[
					    (x < -COUENNE_EPS) ? 0. :  AGGR_MUL);
    else                            return alpha * x + (1. - alpha) * (l + u) / 2.;    // [l,u]
}


/** \brief Constructor. 
 *
 * Get a variable as an argument and set value_ through a call to
 * operator () of that exprAux.
*/

CouenneBranchingObject::CouenneBranchingObject (JnlstPtr jnlst, expression *var, 
						int way, CouNumber brpoint): 
  variable_ (var),
  jnlst_    (jnlst) {

  firstBranch_ =  (way == TWO_LEFT)      ? 0 : 
                 ((way == TWO_RIGHT)     ? 1 : 
                 ((CoinDrand48 () < 0.5) ? 0 : 1));

  CouNumber x = (*variable_) ();//expression::Variable (index_); // current solution

  if (fabs (brpoint) < COUENNE_INFINITY) 
    x = brpoint;

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

  //  assert (fabs (u-l) > COUENNE_EPS);

  CouNumber lb, ub;

  var -> getBounds (lb, ub);

  value_ = midInterval (x, lb, ub);//(*(var -> Lb ())) (), (*(var -> Ub ())) ());
  //  value_ = midInterval (x, expression::Lbound (index_), expression::Ubound (index_));

  //  if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
  jnlst_ -> Printf (J_DETAILED, J_BRANCHING, 
		    "Branch: x%-3d will branch on %g (at %g) [%g,%g]; firstBranch_ = %d\n", 
		    variable_ -> Index (),
		    value_,
		    x, lb, ub,
		    //expression::Variable (index_),
		    //expression::Lbound   (index_),
		    //expression::Ubound   (index_),
		    firstBranch_);
}


/** \brief Execute the actions required to branch, as specified by the
 *	   current state of the branching object, and advance the
 *         object's state.
 *
 *         Returns change in guessed objective on next branch
 */

double CouenneBranchingObject::branch (OsiSolverInterface * solver) {

  // way = 0 if "<=" node, 
  //       1 if ">=" node

  int 
    way   = (!branchIndex_) ? firstBranch_ : !firstBranch_,
    index = variable_ -> Index ();

  bool integer = variable_ -> isInteger ();

  CouNumber
    l    = solver -> getColLower () [index],
    u    = solver -> getColUpper () [index],
    brpt = value_;

  if (jnlst_->ProduceOutput(J_DETAILED, J_BRANCHING)) {
    if (way) {
      if      (value_ < l)             
	jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		       "Nonsense up-br: [ %.8f ,(%.8f)] -> %.8f\n", l,u,value_);
      else if (value_ < l+COUENNE_EPS) 
	jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		       "## WEAK  up-br: [ %.8f ,(%.8f)] -> %.8f\n", l,u,value_);
    } else {
      if      (value_ > u)             
	jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		       "Nonsense dn-br: [(%.8f), %.8f ] -> %.8f\n", l,u,value_);
      else if (value_ > u+COUENNE_EPS) 
	jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		       "## WEAK  dn-br: [(%.8f), %.8f ] -> %.8f\n", l,u,value_);
    }
  }

  if (brpt < l) brpt = l;
  if (brpt > u) brpt = u;

  if (!way) solver -> setColUpper (index, integer ? floor (brpt) : brpt); // down branch
  else      solver -> setColLower (index, integer ? ceil  (brpt) : brpt); // up   branch

  // TODO: apply bound tightening to evaluate change in dual bound

  jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "Branching: x%-3d %c= %g\n", 
  	  index, way ? '>' : '<', integer ? (way ? ceil : floor) (brpt) : brpt);

  branchIndex_++;
  return 0.; // estimated change in objective function
}
