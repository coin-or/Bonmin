/*
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100


/// reoptimize and change bound of a variable if needed
bool updateBound (OsiSolverInterface *csi, /// interface to use as a solver
		  int sense,               /// 1: minimize, -1: maximize
		  CouNumber &bound,        /// bound to be updated
		  bool isint) {            /// is this variable integer

  csi -> setObjSense (sense);
  csi -> setDblParam (OsiDualObjectiveLimit, 
		      (sense == 1) ? COUENNE_INFINITY : -COUENNE_INFINITY); 

  csi -> setDblParam (OsiPrimalObjectiveLimit, bound);
  csi -> resolve ();

  if (csi -> isProvenOptimal ()) {

    double opt = csi -> getObjValue ();

    if (sense > 0) {if (opt > bound + OBBT_EPS) {bound = (isint ? ceil  (opt) : opt); return true;}}
    else           {if (opt < bound - OBBT_EPS) {bound = (isint ? floor (opt) : opt); return true;}}
  }

  return false;
}


/// visit all variables (in one sense) and apply OBBT to each
int obbt_stage (const CouenneCutGenerator *cg,
		OsiSolverInterface *csi,
		OsiCuts &cs,
		t_chg_bounds *chg_bds,
		Bonmin::BabInfo * babInfo,
		int sense) {

  CouenneProblem *p = cg -> Problem ();

  int ncols   = csi -> getNumCols (),
      objind  = p -> Obj (0) -> Body () -> Index (),
      nImprov = 0;

  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  // for all (original+auxiliary) variables x_i,
  ////////////////////////////////////////////////////
  for (int i=0; i<ncols; i++)
    //  for (int i=0; i<problem_ ->nVars(); i++) 

    if ((i != objind) && 
	((sense > 0) && !(chg_bds [i].lower & EXACT) || 
	 (sense < 0) && !(chg_bds [i].upper & EXACT))) { // do not improve objective's bounds

      int nOrig  = p -> nVars ();

      bool isInt = (i < nOrig) ? 
	(p -> Var (i)       -> isInteger ()) :
	(p -> Aux (i-nOrig) -> isInteger ());

      objcoe [i] = 1;
      csi -> setObjective (objcoe);

      CouNumber &bound = (sense == 1) ? (p -> Lb (i)) :	(p -> Ub (i)); 

      // m{in,ax}imize xi on csi

      if (updateBound (csi, sense, bound, isInt)) {

	const double *sol = csi -> getColSolution ();

	if (sense==1) {csi -> setColLower (i, bound); chg_bds [i].lower |= CHANGED | EXACT;}
	else          {csi -> setColUpper (i, bound); chg_bds [i].upper |= CHANGED | EXACT;}

	// check value and bounds of other variables

	for (int j=0; j<ncols; j++) 
	  if ((j!=i) && (j!=objind)) {

	    if (sol [j] <= p -> Lb (j) + COUENNE_EPS) chg_bds [j].lower |= EXACT;
	    if (sol [j] >= p -> Ub (j) - COUENNE_EPS) chg_bds [j].upper |= EXACT;
	  }

	// re-check considering reduced costs (more expensive)

	CouNumber *redcost = NULL;

	// first, compute reduced cost when c = c - e_i, where e_i is
	// a vector with all zero except a one in position i. This
	// serves as a base to compute modified reduced cost below.

	if (0)
	for (int j=0; j<ncols; j++) 
	  if ((j!=i) && (j!=objind)) {

	    // fake a change in the objective function and compute
	    // reduced cost. If resulting vector is all positive
	    // (negative), this solution is also optimal for the
	    // minimization (maximization) of x_j and the
	    // corresponding chg_bds[j] can be set to EXACT.

	    if (sol [j] <= p -> Lb (j) + COUENNE_EPS) chg_bds [j].lower |= EXACT;
	    if (sol [j] >= p -> Ub (j) - COUENNE_EPS) chg_bds [j].upper |= EXACT;
	  }

	// re-apply bound tightening

	if (!(cg -> boundTightening (*csi, cs, chg_bds, babInfo)))
	  return -1; // tell caller this is infeasible

	nImprov++;
      }

      // restore null obj fun
      objcoe [i] = 0;
    }

  return nImprov;
}

/// Optimality based bound tightening
int CouenneCutGenerator::obbt (OsiSolverInterface *csi,
			       OsiCuts &cs,
			       t_chg_bounds *chg_bds,
			       Bonmin::BabInfo * babInfo) const {
  int nImprov;

  csi -> setIntParam (OsiMaxNumIteration, MAX_OBBT_LP_ITERATION);

  // apply all (row+column) cuts to formulation
  csi -> applyCuts (cs);

  nImprov =  obbt_stage (this, csi, cs, chg_bds, babInfo,  1); // minimization
  nImprov += obbt_stage (this, csi, cs, chg_bds, babInfo, -1); // maximization

  return nImprov;
}
