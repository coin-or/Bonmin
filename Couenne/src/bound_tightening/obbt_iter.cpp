/*
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.hpp>
#include <CouenneProblem.hpp>
#include <CouenneSolverInterface.hpp>

#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100

// TODO: seems like Clp doesn't like large bounds and crashes on
// explicit bounds around 1e200 or so. For now simply use fictitious
// bounds around 1e14. Fix.

//#define DEBUG

/// reoptimize and change bound of a variable if needed
bool obbt_updateBound (CouenneSolverInterface *csi, /// interface to use as a solver
		       int sense,               /// 1: minimize, -1: maximize
		       CouNumber &bound,        /// bound to be updated
		       bool isint) {            /// is this variable integer

  csi -> setDblParam (OsiDualObjectiveLimit, COIN_DBL_MAX); 
  csi -> setDblParam (OsiPrimalObjectiveLimit, (sense==1) ? bound : -bound);
  csi -> setObjSense (1); // always minimize, just change the sign of the variable

  ////////////////////////////////////////////////////////////////////////

  csi -> resolve_nobt (); // this is a time-expensive part, be considerate...

  ////////////////////////////////////////////////////////////////////////

  if (csi -> isProvenOptimal ()) {

    double opt = csi -> getObjValue ();

    if (sense > 0) 
         {if (opt          > bound + OBBT_EPS) {bound = (isint ? ceil (opt) : opt); return true;}}
    else {if ((opt = -opt) < bound - OBBT_EPS) {bound = (isint ? floor(opt) : opt); return true;}}
  }

  return false;
}


/// Iteration on one variable

int obbt_iter (const CouenneCutGenerator *cg, 
	       CouenneSolverInterface *csi, 
	       OsiCuts &cs, 
	       t_chg_bounds *chg_bds, 
	       const CoinWarmStart *warmstart, 
	       Bonmin::BabInfo *babInfo,
	       double *objcoe,
	       int sense, 
	       int index) {

  // TODO: do NOT apply OBBT if this is a variable of the form
  // w2=c*w1, as it suffices to multiply result. More in general, do
  // not apply if w2 is a unary monotone function of w1. Even more in
  // general, if w2 is a unary function of w1, apply bound propagation
  // from w1 to w2 and mark it as exact (depending on whether it is
  // non-decreasing or non-increasing

#ifdef DEBUG
  static int iter = 0;
#endif

  std::set <int> deplist;
  int deplistsize;

  bool issimple = false;

  CouenneProblem *p = cg -> Problem ();

  exprVar *var = p -> Var (index);

  int psense  = p   -> Obj (0) -> Sense (),
      objind  = p   -> Obj (0) -> Body () -> Index (),
      ncols   = csi -> getNumCols (),
      nImprov = 0;


  if ((var -> Type  () == AUX) &&
      ((deplistsize = var -> Image () -> DepList (deplist, STOP_AT_AUX, p)) <= 1)) {

    if (!deplistsize) { // funny, the expression is constant...

      CouNumber value = (*(var -> Image ())) ();

      if (csi -> getColLower () [index] < value - COUENNE_EPS) {
	csi -> setColLower (index, value); 
	chg_bds    [index].lower |= CHANGED | EXACT;
      }
      else chg_bds [index].lower |= EXACT;

      if (csi -> getColUpper () [index] > value + COUENNE_EPS) {
	csi -> setColUpper (index, value); 
	chg_bds    [index].lower |= CHANGED | EXACT;
      }
      else chg_bds [index].upper |= EXACT;

      issimple = true;

    } else { // the expression only depends on one variable, meaning
	     // that bound propagation should be sufficient

      int indInd = *(deplist.begin ());

      //      expression *image = var -> Image ();
      // TODO: write code for monotone functions...

      if // independent variable is exactly bounded in both ways
	((chg_bds [indInd].lower & EXACT) && 
	 (chg_bds [indInd].upper & EXACT) ||
	 // or if this expression is of the form w=cx+d, that is, it
	 // depends on one variable only and it is linear
	 (var -> Image () -> Linearity () <= LINEAR)) {

	issimple = true;

	expression *vl, *vu;

	var -> Image () -> getBounds (vl, vu);

	CouNumber 
	  lb = (*vl) (), 
	  ub = (*vu) ();

	delete vl;
	delete vu;

	if (csi -> getColLower () [index] < lb - COUENNE_EPS) {
	  csi -> setColLower (index, lb); 
	  chg_bds      [index].lower |= CHANGED | EXACT;
	} else chg_bds [index].lower |= EXACT;

	if (csi -> getColUpper () [index] > ub + COUENNE_EPS) {
	  csi -> setColUpper (index, ub); 
	  chg_bds      [index].upper |= CHANGED | EXACT;
	} else chg_bds [index].upper |= EXACT;
      }
    }
  }

  // only improve bounds if
  if (!issimple &&
      ((p -> Var (index) -> Type () == VAR) ||             // it is an original variable 
       (p -> Var (index) -> Multiplicity () > 0)) &&       // or its multiplicity is at least 1
      (p -> Lb (index) < p -> Ub (index) - COUENNE_EPS) && // in any case, bounds are not equal

      ((index != objind) || // this is not the objective

       // or it is, so we use it for re-solving
       ((sense ==  1) && (psense == MINIMIZE) && !(chg_bds [index].lower & EXACT)) ||
       ((sense == -1) && (psense == MAXIMIZE) && !(chg_bds [index].upper & EXACT)))) {

    bool isInt = (p -> Var (index) -> isInteger ());

    objcoe [index] = sense;
    csi -> setObjective (objcoe);
    csi -> setObjSense (1); // minimization

    // TODO: Use something else!
#if 0
    for (int iv=0; iv<csi->getNumCols (); iv++) {
      if (fabs (csi -> getColLower () [iv]) > 1e7) csi -> setColLower (iv, -1e14);
      if (fabs (csi -> getColUpper () [iv]) > 1e7) csi -> setColUpper (iv,  1e14);
    }
#endif

    CouNumber &bound = 
      (sense == 1) ? 
      (p -> Lb (index)) : 
      (p -> Ub (index)); 

    // m{in,ax}imize xi on csi

#ifdef DEBUG

    printf ("m%simizing x%d [%g,%g] %c= %g",
	    (sense==1) ? "in" : "ax", index, p -> Lb (index), p -> Ub (index),
	    (sense==1) ? '>'  : '<',  bound); fflush (stdout);

    char fname [20];
    sprintf (fname, "m%s_w%03d_%03d", (sense == 1) ? "in" : "ax", index, iter);
    printf ("\nwriting %s\n", fname);
    csi -> writeLp (fname);
#endif

    csi -> setWarmStart (warmstart);

    if (obbt_updateBound (csi, sense, bound, isInt)) {

      if (p -> bestSol ()) {
	if (sense == 1) {
	  if ((p -> Lb (index) < p -> bestSol () [index]) && 
	      (bound       > COUENNE_EPS + p -> bestSol () [index]))
	    printf ("#### OBBT error on x%d: lb = %g, opt = %g, new lb = %g\n", 
		    index, p -> Lb (index), p -> bestSol () [index], bound);
	} else {
	  if ((p -> Ub (index) > p -> bestSol () [index]) && 
	      (bound       < COUENNE_EPS + p -> bestSol () [index]))
	    printf ("#### OBBT error on x%d: ub = %g, opt = %g, new ub = %g\n", 
		    index, p -> Ub (index), p -> bestSol () [index], bound);
	}
      }

      // more conservative, only change (and set CHANGED) if improve

      if (sense==1)
	if (csi -> getColLower () [index] < bound - COUENNE_EPS) {
#ifdef DEBUG
	  printf ("l_%d: %g --> %g\n", index, csi -> getColLower () [index], bound);
#endif 
	  csi -> setColLower (index, bound); 
	  chg_bds      [index].lower |= CHANGED | EXACT;
	} else chg_bds [index].lower |= EXACT;
      else
	if (csi -> getColUpper () [index] > bound + COUENNE_EPS) {
#ifdef DEBUG
	  printf ("u_%d: %g --> %g\n", index, csi -> getColUpper () [index], bound);
#endif 
	  csi -> setColUpper (index, bound); 
	  chg_bds      [index].upper |= CHANGED | EXACT;
	} else chg_bds [index].upper |= EXACT;

      /*
      if (sense==1) {csi -> setColLower (index, bound); chg_bds [index].lower |= CHANGED | EXACT;}
      else          {csi -> setColUpper (index, bound); chg_bds [index].upper |= CHANGED | EXACT;}
      */

      // check value and bounds of other variables

      const double *sol = csi -> getColSolution ();

      for (int j=0; j<ncols; j++) 
	if ((j!=index) && (j!=objind)) {

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
	if ((j!=index) && (j!=objind)) {

	  // fake a change in the objective function and compute
	  // reduced cost. If resulting vector is all positive
	  // (negative), this solution is also optimal for the
	  // minimization (maximization) of x_j and the corresponding
	  // chg_bds[j].lower (.upper) can be set to EXACT.

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

#ifdef DEBUG
      printf ("XXXXXXXXXXXXXXX OBBT: x_%d: [%g, %g]\n", index, 
	      csi -> getColLower () [index], 
	      csi -> getColUpper () [index]);
#endif

      if (!(cg -> boundTightening (((objind == index) && (sense == psenseI)) ? csi : NULL, 
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

    // TODO: is it, really? And shouldn't you check the opt sense too?
    if ((objind == index) && (csi -> isProvenOptimal ()))
      p -> update (csi -> getColSolution (), NULL, NULL);

    // restore null obj fun
    objcoe [index] = 0;
  }

  return nImprov;
}
