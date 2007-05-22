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


/// Optimality based bound tightening
int CouenneCutGenerator::obbt (const OsiSolverInterface &si,
			       OsiCuts &cs,
			       char *chg_bds,
			       Bonmin::BabInfo * babInfo) const {

  int nImprov = 0, ncols = si.getNumCols ();
  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  // create clone of current solver interface
  OsiSolverInterface *csi = si.clone (true);

  // apply all (row+column) cuts to formulation
  csi -> applyCuts (cs);

  // for all (original+auxiliary) variables x_i,
  for (int i=0; i<ncols; i++) {

    CouNumber opt;
    bool chg = false;

    objcoe [i] = 1;
    csi -> setObjective (objcoe);

    // minimize and then maximize x_i on si.
    csi -> setObjSense  (1); // MINIMIZE
    csi -> resolve ();

    if (csi -> isProvenOptimal () && 
	((opt = csi -> getObjValue ()) > problem_ -> Lb (i) + COUENNE_EPS)) {

      problem_ -> Lb (i) = opt;
      chg = true;
    }

    csi -> setObjSense (-1); // MAXIMIZE
    csi -> resolve ();

    if (csi -> isProvenOptimal () && 
	((opt = csi -> getObjValue ()) < problem_ -> Ub (i) - COUENNE_EPS)) {

      problem_ -> Ub (i) = opt;
      chg = true;
    }

    if (chg) {

      if (!(boundTightening (si, cs, chg_bds, babInfo)))
	return -1; // return value to signal infeasibility

      chg_bds [i] = 1;
      nImprov++;
    }

    // restore null obj fun
    objcoe [i] = 0;
  }

  delete csi;
  return nImprov;
}
