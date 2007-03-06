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


/// Bound tightening

int CouenneCutGenerator::tightenBounds (const OsiSolverInterface &si, 
					char *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // Retrieve, from si, value and bounds of all variables, if not
  // firstcall, otherwise only those of the original ones Update
  // expression structure with x, l, u

  OsiSolverInterface *psi = const_cast <OsiSolverInterface *> (&si);

  CouNumber 
    *xc = const_cast <CouNumber *> (psi -> getColSolution ()),
    *lc = const_cast <CouNumber *> (psi -> getColLower    ()),
    *uc = const_cast <CouNumber *> (psi -> getColUpper    ());

  // update now all variables and bounds

  problem_ -> update (xc, lc, uc);

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. Here,
  // newly enforced branching rules may change dependent auxiliary
  // variables' bounds, in a recursive way. Hence we need to repeat
  // the propagation step as long as at least one bound is modified.

  bool found_one;

  do {

    found_one = false;

    int naux = problem_ -> nAuxs ();

    // check all auxiliary variables for changes in their upper,
    // lower bound, depending on the bound changes of the variables
    // they depend on

    for (register int i = problem_ -> nVars (), j=0; 
	 j < naux; j++) {

      CouNumber ll = (*(problem_ -> Aux (j) -> Lb ())) ();
      CouNumber uu = (*(problem_ -> Aux (j) -> Ub ())) ();

      // check if lower bound got higher    
      if (ll > lc [i+j]) {
	psi -> setColLower (i+j, ll);
	lc [i+j] = ll;
	found_one = true;
      }

      // check if upper bound got lower
      if (uu < uc [i+j]) {
	psi -> setColUpper (i+j, uu);
	uc [i+j] = uu;
	found_one = true;
      }

      if (found_one && chg_bds && !(chg_bds [i+j])) {
	nchg++;
	chg_bds [i+j] = 1;
      }

      expression::update (xc, lc, uc);
    }

  } while (found_one); // repeat as long as at least one bound changed

  // update again 
  problem_ -> update (xc, lc, uc);

  return nchg;
}
