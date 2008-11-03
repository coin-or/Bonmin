/*
 * Name:    CouenneVarObject.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneVarObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneComplObject.hpp"

/// Constructor with information for branching point selection strategy
CouenneVarObject::CouenneVarObject (CouenneProblem *p,
				    exprVar *ref, 
				    Bonmin::BabSetupBase *base, 
				    JnlstPtr jnlst):

  // Do not set variable (expression).
  // This way, no expression-dependent strategy is chosen
  CouenneObject (p, ref, base, jnlst) {

  if (jnlst_ -> ProduceOutput (J_SUMMARY, J_BRANCHING)) {
    printf ("created Variable Object: "); 
    reference_ -> print (); 
    printf (" with %s strategy [clamp=%g, alpha=%g]\n", 
	    (strategy_ == LP_CLAMPED)   ? "lp-clamped" : 
	    (strategy_ == LP_CENTRAL)   ? "lp-central" : 
	    (strategy_ == BALANCED)     ? "balanced"   : 
	    (strategy_ == MIN_AREA)     ? "min-area"   : 
	    (strategy_ == MID_INTERVAL) ? "mid-point"  : 
	    (strategy_ == NO_BRANCH)    ? "no-branching (null infeasibility)" : 
	                                  "no strategy",
	    lp_clamp_, alpha_);
  }
}


/// apply the branching rule 
OsiBranchingObject *CouenneVarObject::createBranch (OsiSolverInterface *si, 
						    const OsiBranchingInformation *info, 
						    int way) const {

  // The infeasibility on an (independent) variable x_i is given by
  // something more elaborate than |w-f(x)|, that is, a function of
  // all infeasibilities of all expressions which depend on x_i.

  problem_ -> domain () -> push 
    (problem_ -> nVars (),
     info -> solution_,
     info -> lower_,
     info -> upper_); // have to alloc+copy

  int bestWay;
  const CouenneObject *criticalObject = NULL; // should create the branchingObject

  CouNumber bestPt = computeBranchingPoint (info, bestWay, criticalObject);

  ///////////////////////////////////////////

  jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, ":::: creating branching on x_%d @%g [%g,%g]\n", 
	  reference_ -> Index (), 
	  info -> solution_ [reference_ -> Index ()],
	  info -> lower_    [reference_ -> Index ()],
	  info -> upper_    [reference_ -> Index ()]);

  OsiBranchingObject *brObj = criticalObject ? 
    criticalObject -> createBranch (si, info, way) :
    new CouenneBranchingObject (si, this, jnlst_, reference_, way, bestPt, doFBBT_, doConvCuts_);

  problem_ -> domain () -> pop ();

  return brObj;
}


/// compute branching point (used in createBranch ())
CouNumber CouenneVarObject::computeBranchingPoint(const OsiBranchingInformation *info,
						  int& bestWay,
						  const CouenneObject *&criticalObject) const {
  criticalObject = NULL;

  if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
    printf ( "---------- computeBRPT for "); 
    reference_ -> print (); 
    printf (" [%g,%g]\n", 
	    info -> lower_ [reference_ -> Index ()],
	    info -> upper_ [reference_ -> Index ()]);
  }

  expression *brVar = NULL; // branching variable

  CouNumber
    brdistDn = 0.,
    brdistUp = 0.,
    bestPt   = 0.,
   *brPts    = NULL, // branching point(s)
   *brDist   = NULL, // distances from LP point to each new convexification
    maxdist = - COIN_DBL_MAX;

  bool chosen = false;

  bestWay = TWO_LEFT;
  int whichWay = TWO_LEFT,
    index = reference_ -> Index ();

  std::set <int> deplist = problem_ -> Dependence () [index];

  for (std::set <int>::iterator i = deplist.begin (); i != deplist.end (); ++i) {

    CouenneObject *obj = problem_ -> Objects () [*i];

    CouNumber improv = 0.;

    assert (obj -> Reference ());

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      printf ("  ** "); 
      obj -> Reference ()             -> print (); 
      if (reference_ -> Image ()) {printf (" := "); obj -> Reference () -> Image () -> print ();}
      printf ("\n");
    }

    if (obj -> Reference ()) {
      if (obj -> Reference () -> Image ())
	improv = obj -> Reference () -> Image ()
	  -> selectBranch (obj, info,                      // input parameters
			   brVar, brPts, brDist, whichWay); // result: who, where, distances, direction
      else {
	brVar = obj -> Reference ();
	brPts  = (double *) realloc (brPts, sizeof (double)); 
	brDist = (double *) realloc (brDist, 2 * sizeof (double)); 

	double point = info -> solution_ [obj -> Reference () -> Index ()];

	*brPts = point;
	improv = 0.;

	if (point > floor (point)) {improv =                  brDist [0] = point - floor (point);}
	if (point < ceil  (point)) {improv = CoinMin (improv, brDist [1] = ceil (point) - point);}

	point -= floor (point);
	whichWay = (point < 0.45) ? TWO_LEFT : (point > 0.55) ? TWO_RIGHT : TWO_RAND;
      }
    }

    if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
      printf ("  --> "); 
      if (brVar) brVar -> print (); 
      printf (" at %g, improv %g <%g>, indices = %d,%d\n", 
	      *brPts, improv, maxdist, index, brVar -> Index ());
    }

    if (brVar &&
	(brVar -> Index () == index) &&    // it's us!
	(fabs (improv) > maxdist) &&       // this branching seems to induce a higher improvement
	(fabs (*brPts) < COU_MAX_COEFF)) { // and branching pt is limited

      criticalObject = problem_ -> Objects () [*i]; // set this object as the branch creator

      brdistDn = brDist [0];
      brdistUp = brDist [1];
      chosen = true;
      bestPt = *brPts;
      maxdist = improv;
      bestWay = whichWay;
    }
  }

  // no hits on this VarObject's variable, that is, this variable was
  // never chosen 

  if (!chosen) {

    bestPt = info -> solution_ [index];
    brVar  = reference_;

    CouNumber 
      l     = info -> lower_ [index], 
      u     = info -> upper_ [index],
      width = lp_clamp_ * (u-l);

    switch (strategy_) {

    case CouenneObject::LP_CLAMPED:
      bestPt = CoinMax (l + width, CoinMin (bestPt, u - width));
      break;
    case CouenneObject::LP_CENTRAL: 
      if ((bestPt < l + width) || (bestPt > u - width))
	bestPt = (l+u)/2;
      break;
    case CouenneObject::MID_INTERVAL: 
    default: 
      // all other cases (balanced, min-area)
      bestPt = midInterval (bestPt, l, u);

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
	if (CoinMin (fabs (bestPt - l), fabs (bestPt - u)) < 1e-3) {
	  printf ("computed failsafe %g [%g,%g] for ", bestPt, l,u);
	  reference_ -> print (); printf ("\n");
	}
      }
      break;
    }

    brPts  = (double *) realloc (brPts, sizeof (double));
    *brPts = bestPt;

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
	printf ("  ::: failsafe:  %g [%g,%g] for ", 
		bestPt, info -> lower_ [index], info -> upper_ [index]); 
	reference_ -> print ();
	printf ("\n");
      }

  } else {

      if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
	if (CoinMin (fabs (bestPt - info -> lower_ [index]), 
		     fabs (bestPt - info -> upper_ [index])) < 1e-3) {
	  printf ("  computed %g [%g,%g] for ", 
		  bestPt, info -> lower_ [index], info -> upper_ [index]); 
	  reference_ -> print ();
	  printf ("\n");
	}
      }
  }

  if (pseudoMultType_ == PROJECTDIST) {

    if (chosen) {
      downEstimate_ = brdistDn;
      upEstimate_   = brdistUp;
    } 
    else downEstimate_ = upEstimate_ = 1.;
  }

  if (brPts)  free (brPts);
  if (brDist) free (brDist);

  return bestPt;
}


#define TOL 0.

/// fix nonlinear coordinates of current integer-nonlinear feasible solution
double CouenneVarObject::feasibleRegion (OsiSolverInterface *solver, 
					 const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  assert (index >= 0);

  double val = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val-TOL);
  solver -> setColUpper (index, val+TOL);

  return 0.;
}


/// are we on the bad or good side of the expression?
bool CouenneVarObject::isCuttable () const {

  const std::set <int>                &deplist = problem_ -> Dependence () [reference_ -> Index ()];
  const std::vector <CouenneObject *> &objects = problem_ -> Objects ();

  for (std::set <int>::const_iterator depvar = deplist. begin ();
       depvar != deplist. end (); ++depvar)
    if (!(objects [*depvar] -> isCuttable ()))
      return false;

  return (!(reference_ -> isInteger ()));
}

