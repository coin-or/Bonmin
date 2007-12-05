/*
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100

///
///

int obbt_iter (const CouenneCutGenerator *cg, 
	       CouenneSolverInterface *csi, 
	       OsiCuts &cs, 
	       t_chg_bounds *chg_bds, 
	       const CoinWarmStart *warmstart, 
	       Bonmin::BabInfo *babInfo,
	       double *objcoe,
	       int sense, 
	       int index);


///
///

int call_iter (const CouenneCutGenerator *cg, 
	       CouenneSolverInterface *csi, 
	       OsiCuts &cs, 
	       t_chg_bounds *chg_bds, 
	       const CoinWarmStart *warmstart, 
	       Bonmin::BabInfo *babInfo,
	       double *objcoe,
	       enum nodeType type,
	       int sense) {

  int ncols   = csi -> getNumCols (),
      nimprov = 0;

  CouenneProblem *p = cg -> Problem ();

  for (int ii=0; ii<ncols; ii++) {
    int i = p -> evalOrder (ii);

    if (p -> Var (i) -> Type () == type) {

      int ni = obbt_iter (cg, csi, cs, chg_bds, warmstart, babInfo, objcoe, sense, i);

      if (ni < 0) return ni;
      nimprov += ni;
    }
  }

  return nimprov;
}


/// Optimality based bound tightening

int CouenneCutGenerator::obbt (CouenneSolverInterface *csi,
			       OsiCuts &cs,
			       t_chg_bounds *chg_bds,
			       Bonmin::BabInfo * babInfo) const {

  // set large bounds to infinity (as per suggestion by JJF)

  int ncols = csi -> getNumCols ();
  const double *lb = csi -> getColLower (),
               *ub = csi -> getColUpper ();

  double inf = csi -> getInfinity ();

  for (int i=ncols; i--;) {
    if (lb [i] < - COUENNE_INFINITY) csi -> setColLower (i, -inf);
    if (ub [i] >   COUENNE_INFINITY) csi -> setColUpper (i,  inf);
  }

  //  csi -> setHintParam (OsiDoDualInResolve, false);

  // setup cloned interface for later use
  csi -> setObjSense (1); // minimization
  csi -> setIntParam (OsiMaxNumIteration, MAX_OBBT_LP_ITERATION);
  csi -> applyCuts (cs);   // apply all (row+column) cuts to formulation
  csi -> initialSolve ();

  const CoinWarmStart *warmstart = csi -> getWarmStart ();

  // improve each bound

  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  csi -> setObjective (objcoe);
  csi -> setObjSense (1);        // minimization

  int nimprov = 0, ni;
 
  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
		  ":::::: OBBT on originals ----------------\n");

  if ((ni = call_iter (this, csi, cs, chg_bds, warmstart, babInfo, objcoe, VAR,  1)) < 0) {
    free (objcoe);
    delete warmstart;
    return ni;
  }

  nimprov += ni;

  if ((ni = call_iter (this, csi, cs, chg_bds, warmstart, babInfo, objcoe, VAR, -1)) < 0) {
    free (objcoe);
    delete warmstart;
    return ni;
  }

  nimprov += ni;

  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
		  ":::::: OBBT on auxiliaries --------------\n");

  if ((ni = call_iter (this, csi, cs, chg_bds, warmstart, babInfo, objcoe, AUX,  1)) < 0) {
    free (objcoe);
    delete warmstart;
    return ni;
  }

  nimprov += ni;

  if ((ni = call_iter (this, csi, cs, chg_bds, warmstart, babInfo, objcoe, AUX, -1)) < 0) {
    free (objcoe);
    delete warmstart;
    return ni;
  }

  Jnlst()->Printf(J_DETAILED, J_BOUNDTIGHTENING,
		  ":::::: ---------------------------------\n");

  free (objcoe);
  delete warmstart;

  return (nimprov + ni);
}
