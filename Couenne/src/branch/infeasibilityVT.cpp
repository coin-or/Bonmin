/*
 * Name:    infeasibilityVT.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Compute violation transfer of a variable
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneVTObject.hpp"


/// compute infeasibility of this variable, |w - f(x)| (where w is the
/// auxiliary variable defined as w = f(x))
double CouenneVTObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  int index = reference_ -> Index ();
  assert (index >= 0);

  if (info -> upper_ [index] - 
      info -> lower_ [index] < COUENNE_EPS)
    return 0.;

  problem_ -> domain () -> push 
    (problem_ -> nVars (),
     info -> solution_, 
     info -> lower_, 
     info -> upper_);

  // get set of variable indices that depend on reference_
  const std::set <int> &dependence = problem_ -> Dependence () [index];

  CouNumber
    xcurr    = info -> solution_ [index],
    retval   = 0.,
    lFeas    = xcurr,
    rFeas    = xcurr,
    leanLeft = 0.;

  if (dependence.size () == 0) { // this is a top level auxiliary,
				 // nowhere an independent

    // that means, for VT, that letting w vary and keeping everything
    // else constant we return the difference between current value
    // and value of function at this point

    if (reference_ -> Image ()) {

      CouNumber
	w  = lFeas = rFeas = (*reference_) (),
	fx = ((*reference_ -> Image ())) ();

      retval = fabs (fx - w); // check if this w=f(x) is used nowhere and is feasible

      if      (lFeas > fx) lFeas = fx;
      else if (rFeas < fx) rFeas = fx;
    }

  } else {

    // this appears as independent in all auxs of the "dependence" list
    // check all non-linear objects containing this variable

    for (std::set <int>::const_iterator i = dependence.begin ();
	 i != dependence.end (); ++i) {

      CouNumber 
	left  = xcurr, 
	right = xcurr;

      const CouenneObject &obj = problem_ -> Objects () [*i];

      // check feasibility of nonlinear object
      if (obj. Reference () && 
	  (fabs ((*(obj. Reference () -> Image ())) () -
		 (*(obj. Reference ())) ())
	   >= CoinMin (COUENNE_EPS, feas_tolerance_))) {

	// measure how far have to go, left and right, to restore
	// feasibility
	obj. Reference () -> Image () -> closestFeasible 
	  (reference_, obj. Reference (), left, right);

	if (left  < lFeas) lFeas = left;
	if (right > rFeas) rFeas = right;
      }

      /*printf ("init retval = %g - %g = %g (diff = %g - %g = %g) ", 
	      rFeas, lFeas, rFeas - lFeas,
	      (*(obj. Reference () -> Image ())) (),
	      (*(obj. Reference ())) (),
	      (*(obj. Reference () -> Image ())) () -
	      (*(obj. Reference ())) ());*/
    }

    retval = rFeas - lFeas;

    if (lFeas < info -> lower_ [index]) lFeas = info -> lower_ [index];
    if (rFeas > info -> upper_ [index]) rFeas = info -> upper_ [index];

    if (rFeas - lFeas > COUENNE_EPS) 
      retval = rFeas - lFeas;
  }

  // if delta is null, return 0, this object is feasible
  if (retval < CoinMin (COUENNE_EPS, feas_tolerance_)) {
    problem_ -> domain () -> pop ();
    return 0.;
  }

  //////////////////////////////////////////////////////////////////////

  // object is infeasible. 
  // check how in the middle of the interval we are

  const CouNumber threshold = .2;

  if (retval > COUENNE_EPS) {

    leanLeft  = (xcurr - lFeas) / retval;

    if      (leanLeft <     threshold) way = 0;
    else if (leanLeft > 1 - threshold) way = 1;
  }

  // done computing delta. Now transfer violation on LP relaxation
  //
  // info -> pi_                         duals
  // info -> solution_                   solution
  // info -> lower_                      lower bound
  // info -> upper_                      upper bound
  // info -> solver () -> getNumRows()   # rows
  // info -> numberColumns_              # cols

  // coefficient in objective is in {0,1} and there is only one
  // variable with coefficient 1
  CouNumber vt_delta =
    (reference_ -> Index () == problem_ -> Obj (0) -> Body () -> Index ()) ? 1. : 0.;

  //   const double * elementByColumn_;
  //   /// Column starts
  //   const CoinBigIndex * columnStart_;
  //   /// Column lengths
  //   const int * columnLength_;
  //   /// Row indices
  //   const int * row_;

  //printf ("------------------ vt_delta [%d] [%g,%g] = %g +\n", index, lFeas, rFeas, vt_delta);

  for (int i=0, n_el = info -> columnLength_ [index]; i < n_el; i++) {

    int indRow = info -> columnStart_ [index] + i;

    vt_delta += 
      fabs (info -> pi_ [info -> row_ [indRow]] * 
	    info -> elementByColumn_  [indRow]);

    /*printf ("+ (pi[%d]=%g) * (el[%d]=%g) [=%g] --> vtd = %g\n",
	    info -> row_ [indRow],
	    info -> pi_ [info -> row_ [indRow]], 
	    indRow,
	    info -> elementByColumn_  [indRow],
	    info -> pi_ [info -> row_ [indRow]] * 
	    info -> elementByColumn_  [indRow],
	    vt_delta);*/
  }

  const CouNumber 
    alpha = 0.7,
    beta  = 0.1;

  /*printf ("return %g * %g + %g * %g + %g * %g --> ", 
	  alpha,        fabs (retval * vt_delta),
	  beta,         retval,
	  1-alpha-beta, leanLeft * (1-leanLeft));*/

  retval = 
    alpha          * fabs (retval*vt_delta) +  // violation transfer itself
    beta           * retval +                  // width of feasibility interval
    (1-alpha-beta) * leanLeft * (1-leanLeft);  // how in the middle of the interval x is

  if ((retval > CoinMin (COUENNE_EPS, feas_tolerance_)) &&
      (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING))) {

    printf ("vt-delta is %-10g [", retval); 

    reference_ -> print (); 
    if (dependence.size () == 0) { // if no list, print image
      printf (" := ");
      reference_ -> Image () -> print ();
    }
    printf ("]\n");
  }

  problem_ -> domain () -> pop ();

  // TODO: add term to prefer branching on larger intervals
  return ((retval < CoinMin (COUENNE_EPS, feas_tolerance_)) ? 
	  0. : retval);
}
