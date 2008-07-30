/*
 * Name:    getDisjunctions.cpp
 * Author:  Pietro Belotti
 * Purpose: get nonlinear (and integer?) disjunctions for the problem
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneDisjCuts.hpp"


/// generate one disjunctive cut from one CGLP
int CouenneDisjCuts::getDisjunctions (std::vector <std::pair <OsiCuts *, OsiCuts *> > &disjs, 
				      OsiSolverInterface &si, 
				      OsiCuts &cs, 
				      const CglTreeInfo &info) const {

  OsiBranchingInformation brInfo (&si, 
				  true,   // bool normalSolver,
				  false); // bool copySolution=false);

  branchingMethod_ -> setupList (&brInfo, true); // initialize

  int 
    ncols  = si.getNumCols (),
    retval = COUENNE_FEASIBLE;

  const double 
    *saveLower = CoinCopyOfArray (si.getColLower (), ncols),
    *saveUpper = CoinCopyOfArray (si.getColUpper (), ncols);

  int nobjs = branchingMethod_ -> numberUnsatisfied ();

  if (nobjs) {

    printf ("%d objects\n", nobjs);

    const int *candidates = branchingMethod_ -> candidates ();
    OsiObject **objs = si. objects ();

    // enumerate first (most infeasible, or most promising) objects of
    // the list

    for (int candInd = 0,      num = 0; 
	 (candInd < nobjs) && (num < numDisjunctions_); 
	 candInd++) {

      CouenneObject *cObj = dynamic_cast <CouenneObject *> (objs [candidates [candInd]]);

      if (!cObj ||                                        // not a nonlinear object
	  (cObj -> checkInfeasibility (&brInfo) == 0.) || // nonlinear, but not violated
	  (cObj -> isCuttable ()))                        // we are on the "bad" (concave) side
	continue;

      // TODO: check if we are on the good or bad side for this object
      // (i.e. implement Cuttable() for all expressions...)

      bool     feasLeft = true,  feasRight = true;
      OsiCuts *leftCuts = NULL, *rightCuts = NULL;

      // left branch ////////////////////////////////////////////////////////////

      OsiBranchingObject *brObj = cObj -> createBranch (&si, &brInfo, 0); // down!

      printf ("cand [%d] = %d: x_%d <>= %g\n", candInd, candidates [candInd],
	      cObj -> Reference () -> Index (), brObj -> value ());

      if (brObj -> branch (&si) >= COUENNE_INFINITY) 
	feasLeft = false;                        // left subtree infeasible
      else leftCuts = getSingleDisjunction (si); // feasible, store new bounds in OsiCuts

      // restore initial bounds
      si.setColLower (saveLower);
      si.setColUpper (saveUpper);

      // right branch ////////////////////////////////////////////////////////////

      brObj = cObj -> createBranch (&si, &brInfo, 1); // up!

      if (brObj -> branch (&si) >= COUENNE_INFINITY) 
	feasRight = false;                        // right subtree infeasible
      else rightCuts = getSingleDisjunction (si); // feasible, store new bounds in OsiCuts

      // restore initial bounds
      si.setColLower (saveLower);
      si.setColUpper (saveUpper);

      // done with generating a disjunction. Is any (or both) side infeasible?

      if (feasLeft && feasRight) {

	// feasible, which variables tightened?

	std::pair <OsiCuts *, OsiCuts *> newpair;
	newpair.first  = leftCuts;
	newpair.second = rightCuts;

	disjs.push_back (newpair);

	num++;

      } else if (feasLeft) {

	// apply leftCuts to si and include them in cs

	retval = COUENNE_TIGHTENED;

      } else if (feasRight) {

	// apply rightCuts to si and include them in cs

	retval = COUENNE_TIGHTENED;

      } else {

	retval = COUENNE_INFEASIBLE;
	break;
      }
    }
  }

  delete [] saveLower;
  delete [] saveUpper;

  // disjunctions (of bounding boxes) generated. Do some extra
  // tightening that may result from tightened sides of the
  // disjunctions

  if (retval != COUENNE_INFEASIBLE) {

    // sanity check: for each disjunction, check if left and right
    // part are feasible with (possibly improved) bounds of si --
    // these may have tightened and now exclude either (or both!)
    // sides of a disjunction

    for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjs.begin ();
	 disjI != disjs.begin (); ++disjI) {

      if   (checkDisjSide (si, disjI -> first)  == COUENNE_INFEASIBLE) { // left  side infeasible?
	if (checkDisjSide (si, disjI -> second) == COUENNE_INFEASIBLE) { // right side infeasible?

	  retval = COUENNE_INFEASIBLE;
	  break;

	} else {
	  // incorporate right side in cs
	  retval = COUENNE_TIGHTENED;
	}
      } else if (checkDisjSide (si, disjI -> second) == COUENNE_INFEASIBLE) { // right side

	// incorporate left side in cs
	retval = COUENNE_TIGHTENED;
      }
    }

    // merge tightened disjunction by intersecting si's bounding box
    // with intersections of smallest boxes, one per disjunction,
    // containing each both sides of the disjunction

    const double
      *lower = si.getColLower (),
      *upper = si.getColUpper ();

    for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjs.begin ();
	 disjI != disjs.begin (); ++disjI) {

      CoinPackedVector lowerChg, upperChg;

      getBoxUnion (si, disjI -> first, disjI -> second, lowerChg, upperChg);

      // apply tightened lower bounds to general problem
      const int    *lindices  = lowerChg.getIndices  ();	
      const double *lvalues   = lowerChg.getElements ();

      for (int i=lowerChg.getNumElements(); i--;) {
	register int    ind = *lindices++;
	register double val = *lvalues++;

	if (val > lower [ind]) {si.setColLower (ind, val); retval = COUENNE_TIGHTENED;}
      }

      // apply tightened upper bounds to general problem
      const int    *uindices  = upperChg.getIndices  ();	
      const double *uvalues   = upperChg.getElements ();

      for (int i=upperChg.getNumElements(); i--;) {
	register int    ind = *uindices++;
	register double val = *uvalues++;

	if (val < upper [ind]) {si.setColUpper (ind, val); retval = COUENNE_TIGHTENED;}
      }
    }   
  }

  return retval;
}
