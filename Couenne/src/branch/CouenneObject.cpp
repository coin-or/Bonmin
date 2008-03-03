/*
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneSolverInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneThreeWayBranchObj.hpp"

#include "exprGroup.hpp"
#include "exprQuad.hpp"
#include "lqelems.hpp"

//#define DEBUG

const CouNumber default_alpha = 0.2;
const CouNumber default_clamp = 0.2;

/// Constructor with information for branching point selection strategy
CouenneObject::CouenneObject (exprVar *ref, Bonmin::BabSetupBase *base,
			      JnlstPtr jnlst):
  reference_      (ref),
  strategy_       (MID_INTERVAL),
  jnlst_          (jnlst),
  alpha_          (default_alpha),
  lp_clamp_       (default_clamp),
  feas_tolerance_ (feas_tolerance_default),
  doFBBT_         (true),
  doConvCuts_     (true) {

  assert (ref -> Type () == AUX);

  if (base) {

    std::string s;

    base -> options() -> GetStringValue ("branch_fbbt",      s,"couenne."); doFBBT_    =(s=="yes");
    base -> options() -> GetStringValue ("branch_conv_cuts", s,"couenne."); doConvCuts_=(s=="yes");

    base -> options () -> GetNumericValue ("feas_tolerance", feas_tolerance_, "couenne.");

    std::string brtype;
    base -> options () -> GetStringValue ("branch_pt_select", brtype, "couenne.");

    if      (brtype == "balanced")  strategy_ = BALANCED;
    else if (brtype == "min-area")  strategy_ = MIN_AREA;
    else if (brtype == "no-branch") strategy_ = NO_BRANCH;
    else if (brtype == "mid-point") {
      strategy_ = MID_INTERVAL;
      base -> options () -> GetNumericValue ("branch_midpoint_alpha", alpha_, "couenne.");
    }

    // accept options for branching rules specific to each operator

    std::string br_operator = "";

    switch (ref -> Image () -> code ()) {

    case COU_EXPRPOW: {

      // begin with default value in case specific exponent are not given
      base -> options () -> GetStringValue ("branch_pt_select_pow", brtype, "couenne.");

      CouNumber expon = ref -> Image () -> ArgList () [1] -> Value ();

      if      (fabs (expon - 2.) < COUENNE_EPS) br_operator = "sqr";
      else if (fabs (expon - 3.) < COUENNE_EPS) br_operator = "cube";
      else if (expon             < 0.)          br_operator = "negpow";
      else                                      br_operator = "pow";
    } break;

    case COU_EXPRMUL: 
      br_operator = (ref -> Image () -> ArgList () [0] -> Index () !=
		     ref -> Image () -> ArgList () [1] -> Index ()) ?
	"prod" : "sqr";
      break;
    case COU_EXPRINV: br_operator = "negpow"; break;
    case COU_EXPRDIV: br_operator = "div"; break;
    case COU_EXPRLOG: br_operator = "log"; break;
    case COU_EXPREXP: br_operator = "exp"; break;
    case COU_EXPRSIN:
    case COU_EXPRCOS: br_operator = "trig"; break;
    default:;
    }

    if (br_operator != "") {
      // read option
      char select [40], sel_clamp [40];
      sprintf (select,    "branch_pt_select_%s", br_operator.c_str ());
      sprintf (sel_clamp, "branch_lp_clamp_%s",  br_operator.c_str ());
      base -> options () -> GetStringValue (select, brtype, "couenne.");
      base -> options () -> GetNumericValue (sel_clamp, lp_clamp_, "couenne.");

      if      (brtype == "balanced")    strategy_ = BALANCED;
      else if (brtype == "lp-clamped")  strategy_ = LP_CLAMPED;
      else if (brtype == "lp-central")  strategy_ = LP_CENTRAL;
      else if (brtype == "min-area")    strategy_ = MIN_AREA;
      else if (brtype == "no-branch")   strategy_ = NO_BRANCH;
      else if (brtype == "mid-point") {
	strategy_ = MID_INTERVAL;
	base -> options () -> GetNumericValue ("branch_midpoint_alpha", alpha_, "couenne.");
      }
    }
  }

  if (jnlst_->ProduceOutput(J_SUMMARY, J_BRANCHING)) {

    printf ("created object: "); reference_ -> print (); 
    printf (" := "); reference_ -> Image () -> print ();
    printf (" with %s strategy [clamp=%g, alpha=%g]\n", 
	    (strategy_ == LP_CLAMPED)   ? "lp-clamped" : 
	    (strategy_ == LP_CENTRAL)   ? "lp-central" : 
	    (strategy_ == BALANCED)     ? "balanced"   : 
	    (strategy_ == MIN_AREA)     ? "min-area"   : 
	    (strategy_ == MID_INTERVAL) ? "mid-point"  : "no-branching (null infeasibility)",
	    lp_clamp_, alpha_);
  }
}


/// Copy constructor
CouenneObject::CouenneObject (const CouenneObject &src):
  reference_  (src.reference_),
  strategy_   (src.strategy_),
  jnlst_      (src.jnlst_),
  alpha_      (src.alpha_),
  lp_clamp_   (src.lp_clamp_),
  feas_tolerance_ (src.feas_tolerance_),
  doFBBT_     (src.doFBBT_),
  doConvCuts_ (src.doConvCuts_) {}


#define TOL 0.

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

    exprGroup::lincoeff &lcoe = e -> lcoeff ();

    for (exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el) {
      int index = el -> first -> Index ();
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

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *si, 
						 const OsiBranchingInformation *info, 
						 int way) const {

  // a nonlinear constraint w = f(x) is violated. The infeasibility is
  // given by something more elaborate than |w-f(x)|, that is, it is
  // the minimum, among the two branching nodes, of the distance from
  // the current optimum (w,x) and the optimum obtained after
  // convexifying the two subproblems. We call selectBranch for the
  // purpose, and save the output parameter into the branching point
  // that should be used later in createBranch.

  CouenneProblem *p = dynamic_cast <CouenneSolverInterface *> (si) -> CutGen () -> Problem ();

  p -> domain () -> push 
    (p -> nVars (),
     info -> solution_,
     info -> lower_,
     info -> upper_);

  CouNumber  *brPts = NULL; // branching point(s)
  expression *brVar = NULL; // branching variable
  int whichWay = 0;

#ifdef DEBUG
  CouNumber improv =  // not used out of debug
#endif

    reference_ -> Image () -> 
    selectBranch (this, info,              // input parameters
		  brVar, brPts, whichWay); // result: who, where, and how to branch

#ifdef DEBUG
  printf ("brpts for "); reference_ -> print (); printf (" := ");
  reference_ -> Image () -> print (); printf (" is on "); brVar -> print ();
  printf (" @ %.12g [%.12g,%.12g]\n", *brPts, 
	  p -> Lb (brVar -> Index ()), 
	  p -> Ub (brVar -> Index ()));
#endif

  // whichWay is IGNORED

#ifdef DEBUG
  if (brVar) {

    if (improv <= COUENNE_EPS) {
      printf ("### warning, infeas = %g for ", improv);
      reference_ -> print (); printf (":=");
      reference_ -> Image () -> print (); printf ("\n");
    }

    int index = brVar -> Index ();
    if (info -> lower_ [index] >= 
	info -> upper_ [index] - COUENNE_EPS) {
      printf ("### warning, tiny bounding box [%g,%g] for x_%d\n", 
	      info -> lower_ [index],
	      info -> upper_ [index], index);
    }
  }
#endif

  // This should make getFixVar() useless if used in exprMul or
  // exprDiv, i.e., the only non-unary operators.

  OsiBranchingObject *brObj = NULL;

  if (brVar) // if applied latest selectBranching

    brObj = new CouenneBranchingObject (jnlst_, brVar, 
					way ? TWO_RIGHT : TWO_LEFT, 
					*brPts, 
					doFBBT_, doConvCuts_);

  else {     // apply default branching rule

    if (jnlst_->ProduceOutput(J_DETAILED, J_BRANCHING)) {
      // we should pipe all output through journalist
      jnlst_->Printf(J_DETAILED, J_BRANCHING, "failsafe branch: ");
      reference_ -> print (std::cout);                              printf (" = ");
      reference_ -> Image () -> print (std::cout); fflush (stdout); printf (" --> branch on ");
      reference_ -> Image () -> getFixVar () -> print (std::cout);  printf ("\n");
    }

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

    CouNumber 
      xr = info -> solution_ [ref_ind],
      lr = info -> lower_    [ref_ind],
      ur = info -> upper_    [ref_ind];

    int index = depvar ? (depvar -> Index ()) : -1;

    if (index >= 0) {

      CouNumber 
	x  = info -> solution_ [index],
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
	brObj = new CouenneBranchingObject (jnlst_, depvar, way ? TWO_RIGHT : TWO_LEFT, x,
					    doFBBT_, doConvCuts_);  
    }

    brObj = new CouenneBranchingObject (jnlst_, reference_, way ? TWO_RIGHT : TWO_LEFT, xr, 
					doFBBT_, doConvCuts_);
  }

  p -> domain () -> pop ();

  if (brPts)
    free (brPts);

  return brObj;
}


/// computes a not-too-bad point where to branch, in the "middle" of an interval
CouNumber CouenneObject::midInterval (CouNumber x, CouNumber l, CouNumber u) const {

  if (u < l + COUENNE_EPS)
    return (0.5 * (l + u));

  if      (x<l) x = l;
  else if (x>u) x = u;

  if   (l < -COUENNE_INFINITY / 10)
    if (u >  COUENNE_INFINITY / 10) return x; // 0.                                    // ]-inf,+inf[
    else                            return ((x < -COUENNE_EPS) ? (AGGR_MUL * (-1+x)) : // ]-inf,u]
					    (x >  COUENNE_EPS) ? 0. : -AGGR_MUL);
  else
    if (u >  COUENNE_INFINITY / 10) return ((x >  COUENNE_EPS) ? (AGGR_MUL *  (1+x)) : // [l,+inf[
					    (x < -COUENNE_EPS) ? 0. :  AGGR_MUL);
    else {                                                                             // [l,u]
      CouNumber point = alpha_ * x + (1. - alpha_) * (l + u) / 2.;
      if      ((point-l) / (u-l) < closeToBounds) point = l + (u-l) * closeToBounds;
      else if ((u-point) / (u-l) < closeToBounds) point = u + (l-u) * closeToBounds;
      return point;
    }
}

/// pick branching point based on current strategy
CouNumber CouenneObject::getBrPoint (funtriplet *ft, CouNumber x0, CouNumber l, CouNumber u) const {

  switch (strategy_) {

  case CouenneObject::MIN_AREA:     return maxHeight   (ft, l, u);
  case CouenneObject::BALANCED:     return minMaxDelta (ft, l, u);
  case CouenneObject::LP_CLAMPED: {
    CouNumber width = lp_clamp_ * (u-l);
    return CoinMax (l + width, CoinMin (x0, u - width));
  }
  case CouenneObject::LP_CENTRAL: {
      CouNumber width = lp_clamp_ * (u-l);
      return ((x0 < l + width) || (x0 > u - width)) ? (l+u)/2 : x0;
  }
  case CouenneObject::MID_INTERVAL: 
  default:                          return midInterval (x0, l, u);
  }
}
