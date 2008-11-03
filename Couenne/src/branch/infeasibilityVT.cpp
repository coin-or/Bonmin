/*
 * Name:    infeasibilityVT.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Compute violation transfer of a variable. See Tawarmalani
 *          and Sahinidis' book or article on MathProg 2002
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

  int indexVar = reference_ -> Index ();
  assert (indexVar >= 0);

  CouNumber tol = CoinMin (COUENNE_EPS, feas_tolerance_);

  // feasible variable
  if (info -> upper_ [indexVar] - 
      info -> lower_ [indexVar] < tol) {
    if (reference_ -> isInteger ()) {
      double point = info -> solution_ [reference_ -> Index ()];
      if (downEstimate_ <       point  - floor (point)) downEstimate_ =       point  - floor (point);
      if (upEstimate_   < ceil (point) -        point)  upEstimate_   = ceil (point) -        point;
      return intInfeasibility (point);
    } else return (upEstimate_ = downEstimate_ = 0.);
  }

  // let Couenne know of current status (variables, bounds)
  problem_ -> domain () -> push 
    (problem_ -> nVars (),
     info -> solution_, 
     info -> lower_, 
     info -> upper_);

  // debug output ////////////////////////////////////////////////////////////////

  if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
    printf ("VT infeas on ");
    reference_ -> print ();
    if (reference_ -> Image ()) { // if no list, print image
      printf (" := ");
      reference_ -> Image () -> print ();
    }
    const std::set <int> &dependence = problem_ -> Dependence () [indexVar];
    if (dependence.size () > 0) {
      printf (" -- ");
      for (std::set <int>::const_iterator i = dependence.begin (); 
	   i != dependence.end (); ++i) {
	problem_ -> Var (*i) -> print ();
	printf (" ");
      }
    }
    printf ("\n");
  }

  // get set of variable indices that depend on reference_
  const std::set <int> &dependence = problem_ -> Dependence () [indexVar];

  CouNumber
    xcurr    = info -> solution_ [indexVar], // current value of variable
    retval   = 0.,                           // will be returned
    lFeas    = xcurr,                        // left-most feasible point with all but x_i constant
    rFeas    = xcurr,                        // right-most ...
    leanLeft = 0.,                           // [0,1],  
    fx       = xcurr,                        // value of expression associated with variable (if aux)
    maxInf   = 0.;                           // max infeasibility over expressions depending on this

  if (reference_ -> Type () == AUX)
    fx = (*(reference_ -> Image ())) ();

  if (dependence.size () == 0) { // this is a top level auxiliary,
				 // nowhere an independent

    // that means, for VT, that letting w vary and keeping everything
    // else constant we return the difference between current value
    // and value of function at this point

    // these variables may also be disabled in BonCouenneSetup.cpp

    assert (reference_ -> Type () == AUX); // otherwise, this is an isolated variable

    // check if this w=f(x) is used nowhere and is feasible
    retval = upEstimate_ = downEstimate_ = maxInf = checkInfeasibility (info);

  } else {

    // this appears as independent in all auxs of the "dependence" list
    // check all non-linear objects containing this variable

    for (std::set <int>::const_iterator i = dependence.begin ();
	 i != dependence.end (); ++i) {

      const CouenneObject *obj = problem_ -> Objects () [*i];
      assert (obj -> Reference ());

      CouNumber 
	left   = xcurr,
	right  = xcurr,
	infeas = obj -> checkInfeasibility (info);

      if (infeas > maxInf)
	maxInf = infeas;

      // check feasibility of nonlinear object
      if (infeas > 0.) {

	// measure how far have to go, left and right, to restore
	// feasibility (see page 582, Tawarmalani and Sahinidis,
	// MathProg A: 99, pp. 563-591)
	//if (obj. Reference () -> Image ()) // not needed! obj only has nonlinear objects
	obj -> Reference () -> Image () -> closestFeasible 
	  (reference_, obj -> Reference (), left, right);

	if (left  < lFeas) lFeas = left;
	if (right > rFeas) rFeas = right;
      }

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) { // debug output
	jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "[%g,%g] --> %g - %g = %g (diff = %g - %g = %g): ", 
			  left, right, rFeas, lFeas, rFeas - lFeas,
			  (obj -> Reference () -> Image ()) ? 
			  (*(obj -> Reference () -> Image ())) () : 0.,  
			  (*(obj -> Reference ())) (),
			  (obj -> Reference () -> Image ()) ? 
			  (*(obj -> Reference () -> Image ())) () - (*(obj -> Reference ())) () : 0.);
	obj ->Reference () -> print (); 
	if (obj -> Reference () -> Image ()) 
	  {printf (" := "); obj -> Reference() -> Image() -> print();}
	printf ("\n");
      }
    }

    if (lFeas < info -> lower_ [indexVar]) lFeas = info -> lower_ [indexVar];
    if (rFeas > info -> upper_ [indexVar]) rFeas = info -> upper_ [indexVar];

    retval = rFeas - lFeas;

    upEstimate_   = rFeas - xcurr;
    downEstimate_ = xcurr - lFeas;

    // need a nonzero value for the estimate. Although maxInf is
    // related to the expressions that depend on this variable, it is
    // still an estimate of how the point will change
    if (upEstimate_   <= tol) upEstimate_   = maxInf;
    if (downEstimate_ <= tol) downEstimate_ = maxInf;
  }

  //////////////////////////////////////////////////////////////////////

  // object is infeasible. 
  // check how in the middle of the interval we are

  const CouNumber threshold = .5; // can be in [0,0.5]

  if (retval > COUENNE_EPS) {

    leanLeft  = (xcurr - lFeas) / retval;

    if      (leanLeft <     threshold) way = 0;
    else if (leanLeft > 1 - threshold) way = 1;
    else way = (CoinDrand48 () < 0.5) ? 0 : 1;
  }

  // done computing delta. Now transfer violation on LP relaxation ///////////////

  // coefficient in objective is in {0,1} and there is only one
  // variable with coefficient 1
  CouNumber vt_delta =
    (indexVar == problem_ -> Obj (0) -> Body () -> Index ()) ? 1. : 0.;

  for (int i=0, n_el = info -> columnLength_ [indexVar]; i < n_el; i++) {

    int indRow = info -> columnStart_ [indexVar] + i;

    vt_delta += 
      fabs (info -> pi_ [info -> row_ [indRow]] * 
    	    info -> elementByColumn_  [indRow]);

    //     vt_delta += 
    //       info -> pi_ [info -> row_ [indRow]] * 
    //       fabs (info -> elementByColumn_  [indRow]);

    jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "+ (pi[%d]=%g) * (el[%d]=%g) [=%g] --> vtd = %g\n",
		      info -> row_ [indRow],
		      info -> pi_ [info -> row_ [indRow]], 
		      indRow,
		      info -> elementByColumn_  [indRow],
		      info -> pi_ [info -> row_ [indRow]] * 
		      info -> elementByColumn_  [indRow],
		      vt_delta);
  }

  // weights for VT itself and width of interval
  const CouNumber 
    alpha = 1.0,
    beta  = 0.0;

  jnlst_ -> Printf (J_MATRIX, J_BRANCHING, "return %g * %g + %g * %g + %g * %g --> ", 
		    alpha, fabs (retval * vt_delta), beta, retval,
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
	for (std::set <int>::const_iterator i = dependence.begin (); 
	     i != dependence.end (); ++i) {
	  problem_ -> Var (*i) -> print ();
	  printf (" ");
	}
      }
      printf ("]\n");
    } else printf ("feasible...\n");
  }

  problem_ -> domain () -> pop ();

  if ((retval < tol) &&
      (maxInf > tol)) {

    // no improvement with this variable, but need to return nonzero
    // if this is an auxiliary whose value is different from relative
    // expression
    retval = maxInf; 
  }

  if (retval < CoinMin (COUENNE_EPS, feas_tolerance_)) 
    retval = 0.;

#define ALMOST_ZERO 1e-8

  assert ((retval < ALMOST_ZERO && upEstimate_ < ALMOST_ZERO && downEstimate_ < ALMOST_ZERO) ||
	  (retval > ALMOST_ZERO && upEstimate_ > ALMOST_ZERO && downEstimate_ > ALMOST_ZERO));

  return (reference_ -> isInteger ()) ? 
    CoinMax (retval, intInfeasibility (info -> solution_ [reference_ -> Index ()])) :
    retval;
}
