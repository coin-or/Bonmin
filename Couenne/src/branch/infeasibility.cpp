/*
 * Name:    infeasibility.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Infeasibility of a non-linear expression
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneThreeWayBranchObj.hpp"

//#define DEBUG

/// return difference between current value
double CouenneObject::infeasibility (const OsiBranchingInformation *info, 
				     int &whichWay) const {

  if (strategy_ == NO_BRANCH)
    return 0.;

  // whichWay should be set to which branch first (for two-way branching?)
  // if selectBranch not called, choose one at random
  whichWay_ = whichWay = TWO_LEFT;
  brVar_ = NULL;

  // infeasibility is always null for linear expressions
  assert ((reference_ -> Image () -> Linearity () > LINEAR) && 
	  (reference_ -> Multiplicity () > 0));

#ifdef DEBUG
  /*printf ("  infeasibility -------------------\n");
  for (int i=0; i<reference_ -> domain () -> current () -> Dimension (); i++)
    printf ("  %4d %20.4g [%20.4g %20.4g] ---> %20.4g [%20.4g %20.4g]\n", i,
	    reference_ -> domain () -> x  (i),
	    reference_ -> domain () -> lb (i),
	    reference_ -> domain () -> ub (i),
	    info -> solution_ [i],
	    info -> lower_    [i],
	    info -> upper_    [i]);*/
#endif

  reference_ -> domain () -> push 
    (reference_ -> domain () -> current () -> Dimension (),
     info -> solution_,
     info -> lower_,
     info -> upper_);

  // if branched-upon variable has a narrow interval, it is not worth
  // to branch on it

  const double expr = (*(reference_ -> Image ())) (), 
               var  = (*reference_) ();

  CouNumber delta = fabs (var - expr);

  /// avoid branching on (relatively and absolutely) small deltas

  if ((delta                           < COUENNE_EPS) || 
      (delta / (1 + fabs (var + expr)) < COUENNE_EPS)) {

#if BR_TEST_LOG >= 0 && defined DEBUG
    if (reference_ -> Image () -> code () == COU_EXPRLOG) {
      printf ("---- found feasible point on curve: ");
      reference_ -> print (); printf (" := ");
      reference_ -> Image () -> print ();
      printf ("\n");
    }
#elif defined DEBUG
    printf ("----|%+g - %+g| = %+g  (delta=%+g) way %d, ind %d. ",  ////[%.2f,%.2f]
	    var, expr, 
	    //	    expression::Lbound (reference_ -> Index ()),
	    //	    expression::Ubound (reference_ -> Index ()),
	    fabs (var - expr), delta, whichWay, reference_ -> Index ());
    reference_             -> print (std::cout); std::cout << " = ";
    reference_ -> Image () -> print (std::cout); printf ("\n");
#endif

    reference_ -> domain () -> pop ();

    return 0.;
  }

  // a nonlinear constraint w = f(x) is violated. The infeasibility
  // is given by something more elaborate than |w-f(x)|, that is, it
  // is the minimum, among the two branching nodes, of the distance
  // from the current optimum (w,x) and the optimum obtained after
  // convexifying the two subproblems. We call selectBranch for the
  // purpose, and save the output parameter into the branching point
  // that should be used later in createBranch.

  // TODO: how to avoid this part of the code when called from
  // CbcModel::setSolution() ?

#if (BR_TEST_LOG >= 0) && BR_TEST_GRAPH
  {
    static bool brStats = true;
    const int NPTS = 300;
    if (brStats) {
      brStats = false;
      double 
	l = info -> lower_ [BR_TEST_LOG],
	u = info -> upper_ [BR_TEST_LOG];
      // draw function on initial (largest) interval
      for (int i=0; i <= NPTS; i++) {
	double x = l + (u-l) * (double) i / (double) NPTS;
	printf ("%10g %10g # brtest\n", x, log (x));
      }
      printf (" # brtest\n#m=0,S=0 # brtest\n\n");
    }
  }
#endif

  CouNumber improv = reference_ -> Image () -> 
    selectBranch (this, info,                // input parameters
		  brVar_, brPts_, whichWay); // result: who, where, and how to branch

  if (brVar_) {

#ifdef DEBUG
    if (improv <= COUENNE_EPS) {
      printf ("### warning, infeas = %g for ", improv);
      reference_ -> print (); printf (":=");
      reference_ -> Image () -> print (); printf ("\n");
    }
#endif

    if (improv > COUENNE_EPS) delta = improv;
    whichWay_ = whichWay;

    // TODO: test this, if commented gives A LOT of problems in nvs24

    int index = brVar_ -> Index ();

    if (info -> lower_ [index] >= 
	info -> upper_ [index] - COUENNE_EPS) {

#ifdef DEBUG
      printf ("### warning, tiny bounding box [%g,%g] for x_%d\n", 
	      info -> lower_ [index],
	      info -> upper_ [index], index);
#endif

      delta = 0.;
    }
  }

  // This should make getFixVar() useless if used in exprMul or
  // exprDiv, i.e., the only non-unary operators.

#ifdef DEBUG
#if BR_TEST_LOG >= 0
  printf ("Inf |%+g - %+g| = %+g  (delta=%+g) way %d, ind %d. [%.10g,%.10g]",
#else
  printf ("Inf |%+g - %+g| = %+g  (delta=%+g) way %d, ind %d. ",  ////[%.2f,%.2f]
#endif
	  var, expr, 
	  fabs (var - expr), delta, whichWay, brVar_ -> Index ()
#if BR_TEST_LOG >= 0
	  , info -> lower_ [BR_TEST_LOG]
	  , info -> upper_ [BR_TEST_LOG]
#endif
	  );

  reference_             -> print (std::cout); std::cout << " = ";
  reference_ -> Image () -> print (std::cout); printf ("\n");
#endif

  reference_ -> domain () -> pop ();

  return fabs (delta);
}
