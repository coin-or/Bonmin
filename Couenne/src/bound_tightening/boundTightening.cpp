/*
 * Name:    boundTightening.cpp
 * Author:  Pietro Belotti
 * Purpose: tighten bounds prior to convexification cuts
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "BonBabInfos.hpp"

// max # bound tightening iterations
#define MAX_BT_ITER 8
#define THRES_IMPROVED 0


// core of the bound tightening procedure

bool CouenneProblem::btCore (t_chg_bounds *chg_bds) const {

  //////////////////////// Bound propagation and implied bounds ////////////////////

  int   ntightened = 0,
      nbwtightened = 0,
      niter = 0;

  bool first = true;

  do {

    // propagate bounds to auxiliary variables

    //    if ((nbwtightened > 0) || (ntightened > 0))
    ntightened = ((nbwtightened > 0) || first) ? 
      tightenBounds (chg_bds) : 0;

    // implied bounds. Call also at the beginning, as some common
    // expression may have non-propagated bounds

    // if last call didn't signal infeasibility
    nbwtightened = ((ntightened > 0) || ((ntightened==0) && first)) ? impliedBounds (chg_bds) : 0;

    if (first)
      first = false;

    if ((ntightened < 0) || (nbwtightened < 0)) {

      // set infeasibility through a cut 1 <= x0 <= -1
      Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
			  "infeasible node at BT\n");
      return false;
    }

    // continue if EITHER procedures gave (positive) results, as
    // expression structure is not a tree.

  } while (((ntightened > 0) || (nbwtightened > 0)) && 
	   (ntightened + nbwtightened > THRES_IMPROVED) &&
	   (niter++ < MAX_BT_ITER));

  // TODO: limit should depend on number of constraints, that is,
  // bound transmission should be documented and the cycle should stop
  // as soon as no new constraint subgraph has benefited from bound
  // transmission. 
  //
  // BT should be more of a graph spanning procedure that moves from
  // one vertex to another when either tightening procedure has given
  // some result. This should save some time especially
  // w.r.t. applying implied bounds to ALL expressions just because
  // one single propagation was found.

  return true;
}


/// procedure to strengthen variable bounds. Return false if problem
/// turns out to be infeasible with given bounds, true otherwise.

bool CouenneProblem::boundTightening (t_chg_bounds *chg_bds, 
				      Bonmin::BabInfo * babInfo) const {

  Jnlst()->Printf (J_DETAILED, J_BOUNDTIGHTENING,
		   "Feasibility-based Bound Tightening\n");

  int objInd = Obj (0) -> Body () -> Index ();

  /////////////////////// MIP bound management /////////////////////////////////

  if ((objInd >= 0) && babInfo && (babInfo -> babPtr ())) {

    CouNumber UB      = babInfo  -> babPtr () -> model (). getObjValue(),
              LB      = babInfo  -> babPtr () -> model (). getBestPossibleObjValue (),
              primal0 = Ub (objInd), 
              dual0   = Lb (objInd);

    // Bonmin assumes minimization. Hence, primal (dual) is an UPPER
    // (LOWER) bound.
    
    if ((UB < COUENNE_INFINITY) && 
	(UB < primal0 - COUENNE_EPS)) { // update primal bound (MIP)

      //      printf ("updating upper %g <-- %g\n", primal0, primal);
      Ub (objInd) = UB;
      chg_bds [objInd].setUpper(t_chg_bounds::CHANGED);
    }

    if ((LB   > - COUENNE_INFINITY) && 
	(LB   > dual0 + COUENNE_EPS)) { // update dual bound
      Lb (objInd) = LB;
      chg_bds [objInd].setLower(t_chg_bounds::CHANGED);
    }
  }

  return btCore (chg_bds);
}


/// procedure to strengthen variable bounds. Return false if problem
/// turns out to be infeasible with given bounds, true otherwise.
int CouenneProblem::redCostBT (const OsiSolverInterface *psi,
			       t_chg_bounds *chg_bds, 
			       Bonmin::BabInfo * babInfo) const {

  if (!babInfo) 
    return 0;

  int nchanges = 0;

  CouNumber
    UB = babInfo -> babPtr () -> model (). getObjValue(),
    LB = babInfo -> babPtr () -> model (). getBestPossibleObjValue ();

  //////////////////////// Reduced cost bound tightening //////////////////////

  if ((LB > -COUENNE_INFINITY) && (UB < COUENNE_INFINITY)) {

    const double 
      *X  = psi -> getColSolution (),
      *RC = psi -> getReducedCost ();

    int ncols = psi -> getNumCols ();

    for (int i=0; i<ncols; i++) {

      CouNumber
	x  = X  [i],
	rc = RC [i],
	dx = Ub (i) - x;

      if ((rc > COUENNE_EPS) && (dx*rc > (UB-LB))) {
	// can improve bound on variable w_i
	Ub (i) = x + (UB-LB) / rc;
	nchanges++;
	chg_bds [i].setUpper(t_chg_bounds::CHANGED);
      }
    }
  }

  return nchanges;
}
