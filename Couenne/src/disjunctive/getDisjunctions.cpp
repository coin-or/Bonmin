/*
 * Name:    getDisjunctions.cpp
 * Author:  Pietro Belotti
 * Purpose: get nonlinear (and integer?) disjunctions for the problem
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneDisjCuts.hpp"


/// generate all disjunctions given current point
int CouenneDisjCuts::getDisjunctions (std::vector <std::pair <OsiCuts *, OsiCuts *> > &disjs, 
				      OsiSolverInterface &si, 
				      OsiCuts &cs, 
				      const CglTreeInfo &info) const {

  OsiBranchingInformation brInfo (&si, 
				  true,   // bool normalSolver,
				  false); // bool copySolution=false);

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS))
    for (int i=0; i<si.getNumCols (); i++)
      printf ("%3d. %8g  [%8g %8g]\n", i, 
	      brInfo. solution_ [i], 
	      brInfo. lower_    [i], 
	      brInfo. upper_    [i]);

  // get list of best disjunction as though we were branching
  branchingMethod_ -> setupList (&brInfo, true); // initialize

  int 
    ncols  = si.getNumCols (),
    retval = COUENNE_FEASIBLE;

  int nobjs = branchingMethod_ -> numberUnsatisfied ();

  if (nobjs) {

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS)) 
      printf ("---   %d disjunctive objects\n", nobjs);

    const int *candidates = branchingMethod_ -> candidates ();
    OsiObject **objs = si. objects ();

    // enumerate first (most infeasible, or most promising) objects of
    // the list

    for (int candInd = 0,      num = 0; 
	 (candInd < nobjs) && (num < numDisjunctions_); 
	 candInd++) {

      CouenneObject *cObj = dynamic_cast <CouenneObject *> (objs [candidates [candInd]]);

      // TODO: consider also integer objects

      if (!cObj ||                                        // IF    not a nonlinear object
	  (cObj -> checkInfeasibility (&brInfo) == 0.) || //    or nonlinear, but not violated
	  (cObj -> isCuttable ()))                        //    or we are on the "good" (convex) side
	continue;                                         // THEN skip

      bool     feasLeft = true,  feasRight = true;
      OsiCuts *leftCuts = NULL, *rightCuts = NULL;

      const double 
	*saveLower = CoinCopyOfArray (si.getColLower (), ncols),
	*saveUpper = CoinCopyOfArray (si.getColUpper (), ncols);

      OsiBranchingObject *brObj = cObj -> createBranch (&si, &brInfo, 0); // down!

      if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	printf ("---   cand [%d] is %d: x_%d <>= %g [%g,%g]\n", 
	      candInd, 
		candidates [candInd],
		dynamic_cast <CouenneBranchingObject *> (brObj) -> variable () -> Index (), 
		brObj -> value (),
		couenneCG_->Problem()->Lb(dynamic_cast<CouenneBranchingObject*>(brObj)
					  ->variable()->Index ()),
		couenneCG_->Problem()->Ub(dynamic_cast<CouenneBranchingObject*>(brObj)
					  ->variable()->Index ()));

      // left branch ////////////////////////////////////////////////////////////

      if (brObj -> branch (&si) >= COUENNE_INFINITY)
	feasLeft = false;                        // left subtree infeasible
      else leftCuts = getSingleDisjunction (si); // feasible, store new bounds in OsiCuts

      // restore initial bounds
      si.setColLower (saveLower);
      si.setColUpper (saveUpper);

      delete brObj;

      // right branch ////////////////////////////////////////////////////////////

      brObj = cObj -> createBranch (&si, &brInfo, 1); // up!

      if (brObj -> branch (&si) >= COUENNE_INFINITY) 
	feasRight = false;                        // right subtree infeasible
      else rightCuts = getSingleDisjunction (si); // feasible, store new bounds in OsiCuts

      delete brObj;

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

	if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	  printf ("---   disj: infeasible right\n");
	applyColCuts (si, leftCuts); // !!! culprit
	retval = COUENNE_TIGHTENED;

      } else if (feasRight) {

	if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	  printf ("---   infeasible left\n");
	applyColCuts (si, rightCuts); // !!! culprit
	retval = COUENNE_TIGHTENED;

      } else {
	  if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	    printf ("---   infeasible!!!\n");
	retval = COUENNE_INFEASIBLE;
	break;
      }

      delete [] saveLower;
      delete [] saveUpper;
    }
  }

  // disjunctions (of bounding boxes) generated. Do some extra
  // tightening that may result from tightened sides of the
  // disjunctions

  if (retval != COUENNE_INFEASIBLE) {

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS))
      printf ("have %d disjunctions\n", disjs.size ());

    // sanity check: for each disjunction, check if left and right
    // part are feasible with (possibly improved) bounds of si --
    // these may have tightened and now exclude either (or both!)
    // sides of a disjunction

    for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjs.begin ();
	 disjI != disjs.end (); ++disjI) {

      if   (checkDisjSide (si, disjI -> first)  == COUENNE_INFEASIBLE) { // left  side infeasible?
	if (checkDisjSide (si, disjI -> second) == COUENNE_INFEASIBLE) { // right side infeasible?

	  retval = COUENNE_INFEASIBLE;
	  break;

	} else {

	  if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	    printf ("---   infeasible left [2]!\n");
	  applyColCuts (si, disjI -> second);
	  retval = COUENNE_TIGHTENED;
	}
      } else if (checkDisjSide (si, disjI -> second) == COUENNE_INFEASIBLE) { // right side

	if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS))
	  printf ("---   infeasible right [2]!\n");
	applyColCuts (si, disjI -> first);
	retval = COUENNE_TIGHTENED;
      }
    }

    if (retval == COUENNE_INFEASIBLE)
      return retval;

    // merge tightened disjunction by intersecting si's bounding box
    // with intersections of smallest boxes, one per disjunction,
    // containing each both sides of the disjunction

    for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjs.begin ();
	 disjI != disjs.end (); ++disjI) {

      CoinPackedVector lowerChg, upperChg;

      // find smallest box containing two disjunctions
      getBoxUnion (si, disjI -> first, disjI -> second, lowerChg, upperChg);

      if ((lowerChg.getNumElements () > 0) || 
	  (upperChg.getNumElements () > 0)) {

	OsiColCut cut;
	cut.setLbs (lowerChg);
	cut.setUbs (upperChg);
	applyColCuts (si, &cut);

	cs.insert (cut);
      }
    }   
  }

  return retval;
}


/// our own applyColCuts
void CouenneDisjCuts::applyColCuts (OsiSolverInterface &si, OsiCuts *cuts) const {

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("applying cuts to SI:\n");
    // apply column cuts to si
    for (int i = cuts -> sizeColCuts (); i--;)
      cuts -> colCutPtr (i) -> print ();
    printf ("--------------------\n");
  }

  // apply column cuts to si
  for (int i = cuts -> sizeColCuts (); i--;)
    applyColCuts (si, cuts -> colCutPtr (i));
}


/// our own applyColCut, single cut
void CouenneDisjCuts::applyColCuts (OsiSolverInterface &si, OsiColCut *cut) const {

  // apply single column cut to si
  if (jnlst_ -> ProduceOutput (J_VECTOR, J_DISJCUTS)) {
    printf ("tightening bounds: "); 
    cut -> print ();
  }

  const CoinPackedVector
    &lbs = cut -> lbs (),
    &ubs = cut -> ubs ();

  const int    *lind = lbs.getIndices  (), *uind = ubs.getIndices  ();
  const double *lval = lbs.getElements (), *uval = ubs.getElements (),
    *oldLower = si.getColLower (), 
    *oldUpper = si.getColUpper ();

  for (int j=lbs.getNumElements (); j--; lval++, lind++)
    if (*lval > oldLower [*lind] + COUENNE_EPS) 
      si.setColLower (*lind, *lval);

  for (int j=ubs.getNumElements (); j--; uval++, uind++)
    if (*uval < oldUpper [*uind] - COUENNE_EPS) 
      si.setColUpper (*uind, *uval);
}
