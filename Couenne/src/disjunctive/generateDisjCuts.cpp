/*
 * Name:    generateDisjCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: separation method for disjunctive cuts
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include "CoinTime.hpp"

#include "CouenneDisjCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

/// generate disjunctive cuts
void CouenneDisjCuts::generateCuts (const OsiSolverInterface &si, 
				    OsiCuts &cs, 
				    const CglTreeInfo info) const {

  if ((depthStopSeparate_ >= 0) &&        // if -1 no limit on depth
      (info.level > depthStopSeparate_))  // check if too deep for adding these cuts
    return;

  double time = CoinCpuTime ();
  OsiSolverInterface *csi = si.clone ();

  // get disjunctions /////////////////////////////
  //
  // start:
  //
  // consider problem Ax <= a_0, x in [l,u]
  //
  // // A and a_0 may be, depending on depth of BB tree and size of the
  // // problem, from the linearization at the root node (in this case
  // // a clone() is needed at the first iteration and branching rules
  // // should be applied explicitly) or the current LP as passed by
  // // si.
  //
  // Get set of disjunctions (in the form of OsiBranchingObjects) in
  // number limited by depth, problem size, and sorted according to
  // branching method (HotInfo if strong branching?)
  //
  // preprocessing /////////////////////////////////
  //
  // for each disjunction (x_i <= or >= x_i^d) {
  //
  //   1) apply left  disj., bound tightening returns x in [l_1,u_1]
  //   2) apply right disj., bound tightening returns x in [l_2,u_2]
  //
  //   3) if both infeasible, bail out
  //   4) if one infeasible only, save column cut, apply tightening to
  //      si and goto start
  // }
  //
  // CGLP //////////////////////////////////////////
  //
  // for each disjunction (x_i <= or >= x_i^d) {
  //
  //   5) get cuts Bx <= b_0 for left  disj.
  //   6) get cuts Cx <= c_0 for right disj.
  //
  //   7) if both infeasible, bail out
  //   8) if one infeasible only, save column cut, apply tightening to
  //      si and continue
  //
  //   9) generate CGLP:
  //         consider subset of Ax <= a_0: 
  //           a) active only, or
  //           b) {root linearization cuts} union {currently active}, or
  //           c) all cuts (!)
  //
  //   10) solve CGLP
  //
  //   11) add corresponding cut, possibly to CGLP as well?
  // }

  numDisjunctions_ = (info.level < depthLevelling_) ? 
    (int) (csi -> numberObjects () * initDisjPercentage_) :
    (int) (csi -> numberObjects () * initDisjPercentage_ / (2 + info.level - depthLevelling_));

  if (numDisjunctions_ < 1) numDisjunctions_ = 1;

  const int 
    exc_infeasible = 1, 
    max_iterations = 4;

  int result;

  try {

    // get disjunctions //////////////////////////////////////////////////////////////

    std::vector <std::pair <OsiCuts *, OsiCuts *> > disjunctions;
    bool start_over;
    int iterations = 0;

    do { // repeat as long as preprocessing shrinks the whole problem

      start_over = false;

      // preprocess
      result = getDisjunctions (disjunctions, *csi, cs, info);

      if      (result == COUENNE_INFEASIBLE) throw exc_infeasible; // fathom node
      else if (result == COUENNE_TIGHTENED)  start_over = true;    // tightened node

      // separation
      for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	   disjI != disjunctions.begin (); ++disjI) {

	result = separateWithDisjunction (*disjI, *csi, cs, info);

	if      (result == COUENNE_INFEASIBLE) throw exc_infeasible; // fathom node
	else if (result == COUENNE_TIGHTENED) {                      // tightened node
	  start_over = true;
	  break;
	}
      }

      if (start_over) {

	for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	     disjI != disjunctions.begin (); ++disjI) {
	  delete disjI -> first;
	  delete disjI -> second;
	}

	disjunctions.erase (disjunctions.begin (), disjunctions.end ());
      }

    } while (start_over && (++iterations < max_iterations));

    // si contains the tightest bounding box. Use it to update
    // CouenneProblem's bounds AND add to cs



    // CGLP //////////////////////////////////////////////////////////////////////////

    // generate one cut for each disjunction

    for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
	 disjI != disjunctions.begin (); ++disjI)

      if (generateDisjCut (*disjI, *csi, cs, info) == COUENNE_INFEASIBLE) // node can be fathomed
	throw exc_infeasible;
  }

  catch (int exception) {

    if (exception == exc_infeasible) {

      // add infeasible column cut 1 <= x_0 <= -1
      OsiColCut *infeascut = new OsiColCut;
      if (infeascut) {
	int i=0;
	double upper = -1., lower = +1.;
	infeascut -> setLbs (1, &i, &lower);
	infeascut -> setUbs (1, &i, &upper);
	cs.insert (infeascut);
	delete infeascut;
      }
    }
  }

  delete csi;

  if (info.level <= 0) 
    nrootcuts_ = cs.sizeRowCuts ();
  ntotalcuts_ += cs.sizeRowCuts ();

  septime_ += (CoinCpuTime () - time);
}
