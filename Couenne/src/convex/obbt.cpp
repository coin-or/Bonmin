/*
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.hpp>
#include <CouenneProblem.hpp>

#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100

//#define DEBUG

/// reoptimize and change bound of a variable if needed
bool updateBound (OsiSolverInterface *csi, /// interface to use as a solver
		  int sense,               /// 1: minimize, -1: maximize
		  CouNumber &bound,        /// bound to be updated
		  bool isint) {            /// is this variable integer

  csi -> setDblParam (OsiDualObjectiveLimit, COIN_DBL_MAX); 
  csi -> setDblParam (OsiPrimalObjectiveLimit, (sense==1) ? bound : -bound);
  csi -> setObjSense (1); // always minimize, just change the sign of the variable

  csi -> resolve ();

  if (csi -> isProvenOptimal ()) {

    double opt = csi -> getObjValue ();

    if (sense > 0) 
         {if (opt          > bound + OBBT_EPS) {bound = (isint ? ceil (opt) : opt); return true;}}
    else {if ((opt = -opt) < bound - OBBT_EPS) {bound = (isint ? floor(opt) : opt); return true;}}
  }

  return false;
}


/// visit all variables (in one sense) and apply OBBT to each
int obbt_stage (const CouenneCutGenerator *cg,
		OsiSolverInterface *csi,
		OsiCuts &cs,
		t_chg_bounds *chg_bds,
		const CoinWarmStart *warmstart,
		Bonmin::BabInfo * babInfo,
		int sense) {

  static int iter = -1;
  iter++;

  CouenneProblem *p = cg -> Problem ();

  int ncols   = csi -> getNumCols (),
      objind  = p -> Obj (0) -> Body () -> Index (),
      nImprov = 0;

  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  csi -> setObjective (objcoe);
  csi -> setObjSense (1); // minimization

  // for all (original+auxiliary) variables x_i,
  ////////////////////////////////////////////////////

  int psense = p -> Obj (0) -> Sense ();

  for (int ii=0; ii<ncols; ii++) {

    int i = p -> evalOrder (ii);

    // TODO: do not apply OBBT if this is a variable of the form w2 =
    // c * w1 as it suffices to multiply result. More in general, do
    // not apply if w2 is a unary monotone function of w1

    // only improve bounds if
    if (((p -> Var (i) -> Type () == VAR) ||       // it is an original variable 
	 (p -> Var (i) -> Multiplicity () > 0)) && // or its multiplicity is at least 1
	
	((i != objind) || // this is not the objective

	 // or it is, so we use it for re-solving
	 ((sense ==  1) && (psense == MINIMIZE) && !(chg_bds [i].lower & EXACT)) ||
	 ((sense == -1) && (psense == MAXIMIZE) && !(chg_bds [i].upper & EXACT)))

	// in both cases, bounds are not equal
	&& (p -> Lb (i) < p -> Ub (i) - COUENNE_EPS)) {

      bool isInt = (p -> Var (i) -> isInteger ());

      //      objcoe [i] = 1;
      objcoe [i] = sense;
      csi -> setObjective (objcoe);
      csi -> setObjSense (1); // minimization

      CouNumber &bound = (sense == 1) ? (p -> Lb (i)) :	(p -> Ub (i)); 

      // m{in,ax}imize xi on csi

#ifdef DEBUG
      printf ("m%simizing x%d [%g,%g] %c= %g",
	      (sense==1) ? "in" : "ax", i, p -> Lb (i), p -> Ub (i),
	      (sense==1) ? '>'  : '<',  bound); fflush (stdout);

      char fname [20];
      sprintf (fname, "m%s_w%03d_%03d", (sense == 1) ? "in" : "ax", i, iter);
      //printf ("\nwriting %s\n", fname);
      //csi -> writeLp (fname);
#endif

      csi -> setWarmStart (warmstart);

      if (updateBound (csi, sense, bound, isInt)) {

	if (p -> bestSol ()) {
	  if (sense == 1) {
	    if ((p -> Lb (i) < p -> bestSol () [i]) && 
		(bound       > COUENNE_EPS + p -> bestSol () [i]))
	      printf ("#### OBBT error on x%d: lb = %g, opt = %g, new lb = %g\n", 
		      i, p -> Lb (i), p -> bestSol () [i], bound);
	  } else {
	    if ((p -> Ub (i) > p -> bestSol () [i]) && 
		(bound       < COUENNE_EPS + p -> bestSol () [i]))
	      printf ("#### OBBT error on x%d: ub = %g, opt = %g, new ub = %g\n", 
		      i, p -> Ub (i), p -> bestSol () [i], bound);
	  }
	}

	if (sense==1) {csi -> setColLower (i, bound); chg_bds [i].lower |= CHANGED | EXACT;}
	else          {csi -> setColUpper (i, bound); chg_bds [i].upper |= CHANGED | EXACT;}

	// check value and bounds of other variables

	const double *sol = csi -> getColSolution ();

	for (int j=0; j<ncols; j++) 
	  if ((j!=i) && (j!=objind)) {

	    if (sol [j] <= p -> Lb (j) + COUENNE_EPS) chg_bds [j].lower |= EXACT;
	    if (sol [j] >= p -> Ub (j) - COUENNE_EPS) chg_bds [j].upper |= EXACT;
	  }

#if 0
	// re-check considering reduced costs (more expensive)

	CouNumber *redcost = NULL;

	// first, compute reduced cost when c = c - e_i, where e_i is
	// a vector with all zero except a one in position i. This
	// serves as a base to compute modified reduced cost below.

	for (int j=0; j<ncols; j++) 
	  if ((j!=i) && (j!=objind)) {

	    // fake a change in the objective function and compute
	    // reduced cost. If resulting vector is all positive
	    // (negative), this solution is also optimal for the
	    // minimization (maximization) of x_j and the
	    // corresponding chg_bds[j] can be set to EXACT.

	    if (!(chg_bds [j].lower & EXACT)) {
	    }

	    if (!(chg_bds [j].upper & EXACT)) {
	    }
	  }
#endif	

	// re-apply bound tightening -- here WITHOUT reduced cost
	// (first argument =NULL is pointer to solverInterface) as csi
	// is not our problem

	int psenseI = (psense == MINIMIZE) ? 1 : -1;

	if (!(cg -> boundTightening (((objind == i) && (sense == psenseI)) ? csi : NULL, 
				     cs, chg_bds, babInfo))) {
#ifdef DEBUG
	  printf ("##### infeasible, bound tightening after OBBT\n");
#endif
	  return -1; // tell caller this is infeasible
	}

	nImprov++;
      }
#ifdef DEBUG
	printf ("\n");
#endif

      // if we solved the problem on the objective function's
      // auxiliary variable (that is, we re-solved the extended
      // problem), it is worth updating the current point (it will be
      // used later to generate new cuts).
      if ((objind == i) && (csi -> isProvenOptimal ()))
	p -> update (const_cast <CouNumber *> (csi -> getColSolution ()), NULL, NULL);

      // restore null obj fun
      objcoe [i] = 0;
    }
  }

  return nImprov;
}


/// Optimality based bound tightening
int CouenneCutGenerator::obbt (OsiSolverInterface *csi,
			       OsiCuts &cs,
			       t_chg_bounds *chg_bds,
			       Bonmin::BabInfo * babInfo) const {

  // set large bounds to infinity (as per suggestion by JJF)

  int ncols = csi -> getNumCols ();
  const double *lb = csi -> getColLower (),
               *ub = csi -> getColUpper ();

  double inf = csi -> getInfinity ();

  while (ncols--) {
    if (lb [ncols] < - COUENNE_INFINITY) csi -> setColLower (ncols, -inf);
    if (ub [ncols] >   COUENNE_INFINITY) csi -> setColUpper (ncols,  inf);
  }

  //  csi -> setHintParam (OsiDoDualInResolve, false);

  // setup cloned interface for later use
  csi -> setObjSense (1); // minimization
  csi -> setIntParam (OsiMaxNumIteration, MAX_OBBT_LP_ITERATION);
  csi -> applyCuts (cs);   // apply all (row+column) cuts to formulation
  csi -> initialSolve ();

  const CoinWarmStart *warmstart = csi -> getWarmStart ();

  // improve each bound

  int nImprovLB = obbt_stage (this, csi, cs, chg_bds, warmstart, babInfo,  1); // lower
  if (nImprovLB < 0) return nImprovLB;

  int nImprovUB = obbt_stage (this, csi, cs, chg_bds, warmstart, babInfo, -1); // upper
  if (nImprovUB < 0) return nImprovUB;

  return (nImprovLB + nImprovUB);
}
