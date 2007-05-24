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

/// reoptimize and change bound of a variable if needed
bool updateBound (OsiSolverInterface *csi, /// interface to use as a solver
		  int sense,               /// 1: minimize, -1: maximize
		  CouNumber &bound,        /// bound to be updated
		  bool isint) {            /// is this variable integer

  csi -> setObjSense (sense);
  csi -> resolve ();

  if (csi -> isProvenOptimal ()) {

    double opt = csi -> getObjValue ();

    if      ((sense>0) && (opt > bound+OBBT_EPS)) bound = (isint ? ceil (opt) : (opt-OBBT_EPS));
    else if ((sense<0) && (opt < bound-OBBT_EPS)) bound = (isint ? floor(opt) : (opt+OBBT_EPS));
    else return false;

    return true;
  }
  else return false;
}


/// Optimality based bound tightening
int CouenneCutGenerator::obbt (const OsiSolverInterface &si,
			       OsiCuts &cs,
			       char *chg_bds,
			       Bonmin::BabInfo * babInfo) const {
  int 
    nImprov = 0, 
    ncols   = si.getNumCols (),
    objind  = problem_ -> Obj (0) -> Body () -> Index ();

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
  for (int i=0; i<ncols; i++) 

    if (i != objind) { // do not improve objective's bounds

      CouNumber opt;
      bool chg = false;
      int nOrig = problem_ -> nVars ();
      bool isInt = (i < nOrig) ? 
	(problem_ -> Var (i)       -> isInteger ()) :
	(problem_ -> Aux (i-nOrig) -> isInteger ());

      objcoe [i] = 1;
      csi -> setObjective (objcoe);

      // minimize and then maximize x_i on si.

      CouNumber l0 = problem_ -> Lb (i);
      CouNumber u0 = problem_ -> Ub (i);

      if (updateBound (csi,  1, problem_ -> Lb (i), isInt)) {
	csi -> setColLower (i, problem_ -> Lb (i));
	chg = true;
      }
      
      if (updateBound (csi, -1, problem_ -> Ub (i), isInt)) {
	csi -> setColUpper (i, problem_ -> Ub (i));
	chg = true;
      }

      if (chg) {
	/*
	printf ("%d: [%g,%g] --> [%g,%g]\n", i,
		l0, u0, 
		problem_ -> Lb (i),
		problem_ -> Ub (i));
	*/
	if (!(boundTightening (si, cs, chg_bds, babInfo)))
	  return -1; // tell caller this is infeasible

	chg_bds [i] = 1;
	nImprov++;
      }

      // restore null obj fun
      objcoe [i] = 0;
    }

  delete csi;
  return nImprov;
}
