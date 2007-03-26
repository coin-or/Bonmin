/*
 * Name:    tightenBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: bound tightening for current linear relaxation
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

// big enough we should stop due to no more bounds tightened
#define MAX_FWTIGHTEN 100 

/// Bound tightening for auxiliary variables

int CouenneProblem::tightenBounds (const OsiSolverInterface &si, 
				   char *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. Here,
  // newly enforced branching rules may change dependent auxiliary
  // variables' bounds, in a recursive way. Hence we need to repeat
  // the propagation step as long as at least one bound is modified.

  // actually, no more than one iteration should be needed, as
  // auxiliary variables only depend on variables with LOWER index

  for (int i=0; i < MAX_FWTIGHTEN; i++) {

    bool found_one = false;

    int naux = nAuxs ();

    // check all auxiliary variables for changes in their upper,
    // lower bound, depending on the bound changes of the variables
    // they depend on

    for (register int i = nVars (), j=0; 
	 j < naux; j++) {

      CouNumber ll = (*(Aux (j) -> Lb ())) ();

      bool chg = false;

      // check if lower bound got higher    
      if ((ll > - COUENNE_INFINITY + 1) && (ll >= lb_ [i+j] + COUENNE_EPS)) {

	//printf ("update lbound %d: %.10f >= %.10f + %.12f\n", i+j, ll, lb_ [i+j], COUENNE_EPS);
	lb_ [i+j] = ll;
	chg = true;
      }

      CouNumber uu = (*(Aux (j) -> Ub ())) ();

      // check if upper bound got lower
      if ((uu < COUENNE_INFINITY - 1) && (uu <= ub_ [i+j] - COUENNE_EPS)) {

	//printf ("update ubound %d: %.10f <= %.10f - %.12f (%.12f)\n", 
	//	i+j, uu, ub_ [i+j], COUENNE_EPS, uu - ub_ [i+j]);
	ub_ [i+j] = uu;
	chg = true;
      }

      if (chg && chg_bds && !(chg_bds [i+j])) {
	nchg++;
	found_one = true;
	chg_bds [i+j] = 1; 
      }

      // useless if assume expression::[lu]b etc already point to
      // problem::[lu]b
      expression::update (NULL, lb_, ub_);
    }

    if (!found_one) break;
  } 

  return nchg;
}
