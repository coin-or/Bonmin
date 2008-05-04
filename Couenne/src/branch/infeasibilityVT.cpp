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

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {

    printf ("VT infeas on ");
    reference_ -> print ();

    if (reference_ -> Image ()) { // if no list, print image
      printf (" := ");
      reference_ -> Image () -> print ();
    }

    const std::set <int> &dependence = problem_ -> Dependence () [reference_ -> Index ()];

    if (dependence.size () > 0) {

      printf (" -- ");

      for (std::set <int>::iterator i = dependence.begin (); 
	   i != dependence.end (); ++i) {
	problem_ -> Var (*i) -> print ();
	printf (" ");
      }
    }
    printf ("\n");
  }

  int index = reference_ -> Index ();
  assert (index >= 0);

  /*if (info -> upper_ [index] - 
    info -> lower_ [index] < CoinMin (COUENNE_EPS, feas_tolerance_))
    return (upEstimate_ = downEstimate_ = 0.);*/

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
    leanLeft = 0.,
    w        = (*reference_) (),
    fx       = w;

  if (reference_ -> Type () == AUX)
    fx = (*(reference_ -> Image ())) ();

  if (dependence.size () == 0) { // this is a top level auxiliary,
				 // nowhere an independent

    // that means, for VT, that letting w vary and keeping everything
    // else constant we return the difference between current value
    // and value of function at this point

    if (reference_ -> Type () == AUX) {

      lFeas = rFeas = w;

      //retval = fabs (fx - w); // check if this w=f(x) is used nowhere and is feasible

      if      (lFeas > fx) lFeas = fx;
      else if (rFeas < fx) rFeas = fx;
    } // otherwise, this is an isolated variable

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

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {

	jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "[%g,%g] --> %g - %g = %g (diff = %g - %g = %g): ", 
			  left, right, rFeas, lFeas, rFeas - lFeas,
			  (*(obj. Reference () -> Image ())) (), (*(obj. Reference ())) (),
			  (*(obj. Reference () -> Image ())) () - (*(obj. Reference ())) ());

	obj.Reference () -> print (); printf (" := ");
	obj.Reference () -> Image () -> print (); printf ("\n");
      }
    }

    //retval = rFeas - lFeas;

    if (lFeas < info -> lower_ [index]) lFeas = info -> lower_ [index];
    if (rFeas > info -> upper_ [index]) rFeas = info -> upper_ [index];

    //if (rFeas - lFeas > COUENNE_EPS) 
    retval = rFeas - lFeas;
  }

  // if delta is null, return 0, this object is feasible
  //if (retval < CoinMin (COUENNE_EPS, feas_tolerance_)) {

  //jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "2nd step, infeas = 0\n");
  //problem_ -> domain () -> pop ();
  //return (upEstimate_ = downEstimate_ = 0.);
  //}

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

    jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "+ (pi[%d]=%g) * (el[%d]=%g) [=%g] --> vtd = %g\n",
		      info -> row_ [indRow],
		      info -> pi_ [info -> row_ [indRow]], 
		      indRow,
		      info -> elementByColumn_  [indRow],
		      info -> pi_ [info -> row_ [indRow]] * 
		      info -> elementByColumn_  [indRow],
		      vt_delta);
  }

  const CouNumber 
    alpha = 1.0,
    beta  = 0.0;

  jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "return %g * %g + %g * %g + %g * %g --> ", 
		    alpha,        fabs (retval * vt_delta),
		    beta,         retval,
		    1-alpha-beta, leanLeft * (1-leanLeft));

  retval = 
    alpha          * fabs (retval*vt_delta) +  // violation transfer itself
    beta           * retval +                  // width of feasibility interval
    (1-alpha-beta) * leanLeft * (1-leanLeft);  // how in the middle of the interval x is

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {

    if (retval > CoinMin (COUENNE_EPS, feas_tolerance_)) {

      printf ("vt-delta is %-10g [", retval); 

      reference_ -> print (); 

      if (reference_ -> Image ()) { // if no list, print image
	printf (" := ");
	reference_ -> Image () -> print ();
      }

      if (dependence.size () > 0) {

	printf (" -- ");

	for (std::set <int>::iterator i = dependence.begin (); 
	     i != dependence.end (); ++i) {
	  problem_ -> Var (*i) -> print ();
	  printf (" ");
	}
      }

      printf ("]\n");
    } else {
      printf ("feasible...\n");
    }
  }

  if (retval < CoinMin (COUENNE_EPS, feas_tolerance_))
    retval = 0.;

  //upEstimate_ = downEstimate_ = retval;

  problem_ -> domain () -> pop ();

  if (retval < CoinMin (COUENNE_EPS, feas_tolerance_)) {

    // apparently no improvement with this variable, but need to
    // return nonzero if this is an auxiliary whose value is different from 

    if (fabs (fx - w) > CoinMin (COUENNE_EPS, feas_tolerance_))
      retval = .5; // ensure this variable will have lesser priority
		   // than those with nonzero VT infeasibility --
		   // these will be added 1, see below
  } else retval += 1.;

  upEstimate_ = downEstimate_ = retval;
  return retval;
}
