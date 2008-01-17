/*
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneThreeWayBranchObj.hpp"

#include "exprGroup.hpp"
#include "exprQuad.hpp"
#include "lqelems.hpp"


/// Constructor with information for branching point selection strategy
CouenneObject::CouenneObject (exprVar *ref, Bonmin::BabSetupBase *base,
			      JnlstPtr jnlst):

  reference_ (ref),
  brVar_     (NULL), 
  brPts_     (NULL),
  whichWay_  (BRANCH_NONE),
  strategy_  (MID_INTERVAL),
  jnlst_     (jnlst) {

  if (ref -> Type () == VAR) {
    printf ("Couenne error: CouenneObject cannot be defined on original variables\n");
    exit (-1);
  }

  if (base) {

    std::string brtype;
    base -> options () -> GetStringValue ("branch_pt_select", brtype, "couenne.");

    if      (brtype == "balanced")  strategy_ = BALANCED;
    else if (brtype == "min-area")  strategy_ = MIN_AREA;
    else if (brtype == "mid-point") strategy_ = MID_INTERVAL;
  }
}


/// Copy constructor
CouenneObject::CouenneObject (const CouenneObject &src):
  reference_ (src.reference_),
  brVar_     (src.brVar_),
  brPts_     (NULL),
  whichWay_  (src.whichWay_),
  strategy_  (src.strategy_),
  jnlst_     (src.jnlst_) {

  if (src.brPts_) {

    int nbrpts = 
      ((whichWay_ == TWO_RAND) ||
       (whichWay_ == TWO_LEFT) ||
       (whichWay_ == TWO_RIGHT)) ? 1 : 2;

    brPts_ = (CouNumber *) malloc (nbrpts * sizeof (CouNumber));

    while (nbrpts--) 
      brPts_ [nbrpts] = src.brPts_ [nbrpts];
  }
}


#define TOL 0

/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  assert (index >= 0);

  double val = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val-TOL);
  solver -> setColUpper (index, val+TOL);

  expression *expr = reference_ -> Image ();

  // fix all variables upon which this auxiliary depends

  // expr is surely nonlinear, so it's useless to check if it is an
  // exprAux, w1:=w0

  if (expr -> Type () == UNARY) { // unary function

    index = expr -> Argument () -> Index ();

    if (index >= 0) {
      val = info -> solution_ [index];
      solver -> setColLower (index, val-TOL);
      solver -> setColUpper (index, val+TOL);
    }
  }
  else // n-ary function

    if (expr -> Type () == N_ARY) {

      expression ** args = expr -> ArgList ();
      int nargs = expr -> nArgs ();

      for (register int i=0; i < nargs; i++) {

	if ((index = args [i] -> Index()) >= 0) {
	  val = info -> solution_ [index];
	  solver -> setColLower (index, val-TOL);
	  solver -> setColUpper (index, val+TOL);
	}
      }
    }

  // last cases: exprGroup/Quad, must handle the linear/quadratic terms
  if ((expr -> code () == COU_EXPRGROUP) ||
      (expr -> code () == COU_EXPRQUAD)) {

    exprGroup *e = dynamic_cast <exprGroup *> (expr);

    const exprGroup::lincoeff &lcoe = e -> lcoeff ();

    for (int n = lcoe.size(), i=0; n--; i++) {
      //exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el) {
      int index = lcoe [i]. first -> Index ();
      val = info -> solution_ [index];
      solver -> setColLower (index, val-TOL);
      solver -> setColUpper (index, val+TOL);
    }

    // take care of quadratic terms
    if (expr -> code () == COU_EXPRQUAD) {

      exprQuad *e = dynamic_cast <exprQuad *> (expr);

      exprQuad::sparseQ q = e -> getQ ();

      for (exprQuad::sparseQ::iterator row = q.begin (); 
	   row != q.end (); ++row) {

	int xind = row -> first -> Index ();

	val = info -> solution_ [xind];
	solver -> setColLower (xind, val-TOL);
	solver -> setColUpper (xind, val+TOL);

	for (exprQuad::sparseQcol::iterator col = row -> second.begin ();
	     col != row -> second.end (); ++col) {

	  int yind = col -> first -> Index ();

	  val = info -> solution_ [yind];
	  solver -> setColLower (yind, val-TOL);
	  solver -> setColUpper (yind, val+TOL);
	}
      }
    }
  }

  // TODO: better value through one run of btCore ()

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *si, 
						 const OsiBranchingInformation *info, 
						 int way) const {

  bool isint = (brVar_) && brVar_ -> isInteger (); 
  // (brVarInd_ >= 0) && (si -> isInteger (brVarInd_));

  // way has suggestion from CouenneObject::infeasibility(), but not
  // as set in infeasibility, so we use the one stored in member
  // whichWay_
  // AW: We can't do that, way must be used!

  // way = whichWay_;

#if (BR_TEST_LOG >= 0) && BR_TEST_GRAPH

  //#define OPT_X0 1.7873085033688871
    {
      static bool first = true;

      if (first) {
	first = false;
      }

      double 
	l = info -> lower_ [BR_TEST_LOG],
	u = info -> upper_ [BR_TEST_LOG];
      if (//(l <= OPT_X0) && (u >= OPT_X0) &&
	  (reference_ -> Image () -> code () == COU_EXPRLOG))
	printf ("#m=1,S=0 # brtest\n\
%10g %10g # brtest\n\
%10g %10g # brtest\n\
%10g %10g # brtest\n\
%10g %10g # brtest\n\
%10g %10g # brtest\n\
 # brtest\n\
#m=-1,S=0 # brtest\n\
 # brtest\n",
		l, log (l), 
		//		*brPts_, log (*brPts_), 
		*brPts_, log (*brPts_)-CoinMin(u-*brPts_, *brPts_-l)/2, 
		*brPts_, log (*brPts_), 
		*brPts_, log (*brPts_)-CoinMin(u-*brPts_, *brPts_-l)/2, 
		u, log (u));
    }
#endif


  if (brVar_) // if applied latest selectBranching

    switch (way) {

    case TWO_LEFT:
    case TWO_RIGHT:
    case TWO_RAND:
      jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		     "2way Branch x%d at %g [%d] (%d)\n", 
		     brVar_ -> Index (), *brPts_, way, isint);
      return new CouenneBranchingObject (jnlst_, brVar_, way, *brPts_);
    case THREE_LEFT:
    case THREE_CENTER:
    case THREE_RIGHT:
    case THREE_RAND:
      jnlst_->Printf(J_DETAILED, J_BRANCHING, 
		     "3Way Branch x%d @ %g ][ %g [%d] (%d)\n", 
		     brVar_ -> Index (), *brPts_, brPts_ [1], way, isint);
      return new CouenneThreeWayBranchObj (jnlst_, brVar_, brPts_ [0], brPts_ [1], way);
    default: 
      printf ("CouenneObject::createBranch(): way=%d has no sense\n", way);
      exit (-1);
    }

  // if selectBranch returned -1, apply default branching rule

  if (jnlst_->ProduceOutput(J_DETAILED, J_BRANCHING)) {
    // we should pipe all output through journalist
    jnlst_->Printf(J_DETAILED, J_BRANCHING, "CO::createBranch: ");
    reference_ -> print (std::cout);                              printf (" = ");
    reference_ -> Image () -> print (std::cout); fflush (stdout); printf (" --> branch on ");
    reference_ -> Image () -> getFixVar () -> print (std::cout);  printf ("\n");
  }

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

  expression *depvar = reference_ -> Image () -> getFixVar ();

  // Create a two-way branching object according to finiteness of the
  // intervals. For now only check if argument bounds are finite.

  int ref_ind = reference_ -> Index ();

  CouNumber xr = info -> solution_ [ref_ind],
            lr = info -> lower_    [ref_ind],
            ur = info -> upper_    [ref_ind];

  int index = depvar ? (depvar -> Index ()) : -1;

  if (index >= 0) {

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
      return new CouenneBranchingObject (jnlst_, depvar, way, x);  
  }

  return new CouenneBranchingObject (jnlst_, reference_, way, xr);
}
