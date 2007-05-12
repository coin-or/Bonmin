/*
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <stdlib.h>

#include <CouenneObject.hpp>
#include <exprGroup.h>
#include <CouenneBranchingObject.hpp>
#include <CouenneThreeWayBranchObj.hpp>

//#define WEI_INF   1.
//#define WEI_RANK  0.
//#define WEI_MULT  0.


/// return difference between current value
double CouenneObject::infeasibility (const OsiBranchingInformation *info, 
				     int &whichWay) const {

  // whichWay should be set to which branch first (for two-way branching?)

  // infeasibility is always null for linear expressions
  if (reference_ -> Image () -> Linearity () <= LINEAR)
    return 0.;

  expression::update (const_cast <CouNumber *> (info -> solution_),
		      const_cast <CouNumber *> (info -> lower_),
		      const_cast <CouNumber *> (info -> upper_));

  expression *fixvar = reference_ -> Image () -> getFixVar ();
  int index = fixvar -> Index ();

  if (index < 0)
    return 0.;

  // if branched-upon variable has a narrow interval, it is not worth
  // to branch on it

  const double & expr = (*(reference_ -> Image ())) (), 
               & var  = (*reference_) ();

  //expression::Variable (reference_ -> Index ());
  /* if (0) {

    printf ("Inf: = |%.6e - %.6e| = %.6e  ",  ////[%.2f,%.2f]
	    var, expr, 
	    //	    expression::Lbound (reference_ -> Index ()),
	    //	    expression::Ubound (reference_ -> Index ()),
	    fabs (var - expr));
    reference_             -> print (std::cout); std::cout << " = ";
    reference_ -> Image () -> print (std::cout); printf ("\n");
  }  */

  CouNumber delta = fabs (var - expr);

  //  CouNumber l  = info -> lower_ [index],
  //            u  = info -> upper_ [index];
  /*  if (0)
    if ((delta > COUENNE_EPS) &&
	((fabs (u-l) < COUENNE_EPS) ||
	((mymin (fabs (l), fabs (u)) > COUENNE_EPS) && 
	 (fabs (u-l) / mymax (fabs (l), fabs (u)) < COUENNE_EPS)))) {
      //      ((mymin (fabs (lr), fabs (ur)) > COUENNE_EPS) && 
      //       (fabs (ur-lr) / mymax (fabs (lr), fabs (ur)) < COUENNE_EPS)))

      printf (". Inf: = |%.4f - %.4f| = %.4e. w [%.3f,%.3f], x [%.3f,%.3f] = %.4e ",  ////[%.2f,%.2f]
	      var, expr, fabs (var - expr), 
	      info -> lower_ [reference_ -> Index ()],
	      info -> upper_ [reference_ -> Index ()],
	      l, u, u-l);
      reference_             -> print (std::cout); std::cout << " = ";
      reference_ -> Image () -> print (std::cout);
      printf ("\n");
    }
  */
  //printf (" delta=%.9f,l=%.9f,u=%.9f ", delta, l, u);

  /// avoid branching on (relatively) small deltas
  if (delta < COUENNE_EPS)
    /*||
      (fabs (u-l) < COUENNE_EPS) ||
      ((mymin (fabs (l), fabs (u)) > COUENNE_EPS) && 
      (fabs (u-l) / mymax (fabs (l), fabs (u)) < COUENNE_EPS)))*/
    //      ((mymin (fabs (lr), fabs (ur)) > COUENNE_EPS) && 
    //       (fabs (ur-lr) / mymax (fabs (lr), fabs (ur)) < COUENNE_EPS)))
    delta = 0.;

  else {

    // a nonlinear constraint w = f(x) is violated. The infeasibility
    // is given by something more elaborate than |w-f(x)|, that is, it
    // is the minimum, among the two branching nodes, of the distance
    // from the current optimum (w,x) and the optimum obtained after
    // convexifying the two subproblems. We call selectBranch for the
    // purpose, and save the output parameter into the branching point
    // that should be used later in createBranch.

    CouNumber improv = reference_ -> Image () -> selectBranch 
      (reference_, info,             // input parameters
       brVarInd_, brPts_, whichWay); // result: who, where, and how to branch

    if (brVarInd_ >= 0)
      delta = improv;

    // This should make getFixVar() useless if used in exprMul or
    // exprDiv, i.e., the only non-unary operators.

    // TODO: how to avoid this part of the code when called from
    // CbcModel::setSolution() ?

    // make delta a function of the variable's rank and multiplicity

    /*delta =   WEI_INF  * (1. - exp (-delta))
            + WEI_RANK / (1. + fixvar -> rank ())
            + WEI_MULT * (1. - 1. / fixvar -> Multiplicity ());*/
  }

  return delta;
}


/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  if (index < 0) {   // should never happen...
    printf ("Warning, CouenneObject::feasibleRegion: reference_'s index negative\n");
    return 0;
  }

  double val = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val);
  solver -> setColUpper (index, val);

  expression * expr = reference_ -> Image ();

  // fix all variables upon which this auxiliary depends

  if (expr -> Argument ()){ // unary function

    index = expr -> Argument () -> Index ();

    if (index >= 0) {
      val = info -> solution_ [index];
      solver -> setColLower (index, val);
      solver -> setColUpper (index, val);
    }
  }
  else // n-ary function
    if (expr -> ArgList ()) {

      expression ** args = expr -> ArgList ();
      int nargs = expr -> nArgs ();

      for (register int i = 0 ; i < nargs ; i++) {

	if ((index = args [i] -> Index()) >= 0) {
	  val = info -> solution_ [index];
	  solver -> setColLower (index, val);
	  solver -> setColUpper (index, val);
	}
      }
    }

  // last case: exprGroup, the linear terms are not handled
  if (expr -> code () == COU_EXPRGROUP) {

    exprGroup *e = dynamic_cast <exprGroup *> (expr);
    int *indices = e -> getIndices ();

    for (; *indices >= 0; indices++) {
      val = info -> solution_ [*indices];
      solver -> setColLower (*indices, val);
      solver -> setColUpper (*indices, val);
    }
  }

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *si, 
						 const OsiBranchingInformation *info, 
						 int way) const {

  // way has suggestion from CouenneObject::infeasibility()

  bool isint = reference_ -> isInteger ();

  if (brVarInd_ >= 0)
    switch (way) {
    case TWO_LEFT:
    case TWO_RIGHT:
    case TWO_RAND:
      return new CouenneBranchingObject (brVarInd_, way, *brPts_, isint);
    case THREE_LEFT:
    case THREE_CENTER:
    case THREE_RIGHT:
    case THREE_RAND:
      return new CouenneThreeWayBranchObj (brVarInd_, brPts_ [0], brPts_ [1], way, isint);
    default: 
      printf ("CouenneObject::createBranch(): way has no sense\n");
      exit (-1);
    }

  /*if (0) {
    printf ("CO::createBranch: ");
    reference_ -> print (std::cout);
    printf (" = ");
    reference_ -> Image () -> print (std::cout);
    printf (" --> branch on ");
    reference_ -> Image () -> getFixVar () -> print (std::cout);
    printf ("\n");
    }*/

  // constructor uses actual values of variables and bounds, update them
  expression::update (const_cast <CouNumber *> (info -> solution_),
		      const_cast <CouNumber *> (info -> lower_),
		      const_cast <CouNumber *> (info -> upper_));

  // change the value of delta to reflect the branching operations
  // that will take place. This implies repeatedly faking generation
  // of convexification cuts for different branching points until we
  // have a good branching point. 
  //
  // The infeasibility returned is the minimum of the distances from
  // the current point to the two new convexifications, which is the
  // function that we want to maximize.

  //  CouNumber delta = reference_ -> Image () -> BranchGain (reference_, info);

  // operator-dependent branching object
  //  OsiBranchingObject *opobj = reference_ -> Image () -> BranchObject (reference_, info);

  //  if (opobj)
  //    return opobj;

  expression *depvar = reference_ -> Image () -> getFixVar ();
  int index;

  // Create a two- or three-way branching object according to
  // finiteness of the intervals. For now only check if argument
  // bounds are finite.

  // TODO: check if function is finite within the interval (if so,
  // two-way, if not, three-way). Could be operator-dependent, that
  // is,
  //
  // return reference_ -> BranchObject (brVar_);

  int ref_ind = reference_ -> Index ();

  CouNumber xr = info -> solution_ [ref_ind],
            lr = info -> lower_    [ref_ind],
            ur = info -> upper_    [ref_ind];

  if (depvar && ((index = depvar -> Index ()) >= 0)) {

    CouNumber x  = info -> solution_ [index],
              l  = info -> lower_    [index],
              u  = info -> upper_    [index];
    /*
    if (((x-l > COUENNE_LARGE_INTERVAL) &&
	 (u-x > COUENNE_LARGE_INTERVAL)) 
	|| 
	(((x-l > COUENNE_LARGE_INTERVAL) ||
	  (u-x > COUENNE_LARGE_INTERVAL)) && 
	 ((x-l < COUENNE_NEAR_BOUND) ||
	  (u-x < COUENNE_NEAR_BOUND))))
      return new CouenneThreeWayBranchObj (depvar, x, l, u);
    */

    if (((fabs (x-l) > COUENNE_EPS) &&
	 (fabs (u-x) > COUENNE_EPS) &&
	 (fabs (u-l) > COUENNE_EPS))
	|| (fabs (xr-lr) < COUENNE_EPS)
	|| (fabs (ur-xr) < COUENNE_EPS)
	|| (fabs (ur-lr) < COUENNE_EPS))
      return new CouenneBranchingObject (index, TWO_RAND, x, depvar -> isInteger ());
  }

  return new CouenneBranchingObject (ref_ind, TWO_RAND, xr, reference_ -> isInteger ());
}
