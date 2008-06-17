/*
 * Name:    CouenneLPtightenBounds.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tighten LP bounds on all variables (including continuous)
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

// Tighten bounds - lightweight. Returns -1 if infeasible, otherwise
// number of variables tightened.
int CouenneSolverInterface::tightenBounds (int lightweight) {

  // not yet...
  return OsiClpSolverInterface::tightenBounds (lightweight);

  int 
    ncols = getNumCols (),
    nTightened;

  double 
    *oldLower = new double [ncols],
    *oldUpper = new double [ncols];

  CoinCopyN (getColLower (), ncols, oldLower);
  CoinCopyN (getColUpper (), ncols, oldUpper);

  nTightened = OsiClpSolverInterface::tightenBounds (lightweight);

  if (nTightened > 0) {

    // something was tightened. Run an extra btCore "por si las
    // moscas" (just in case)

    const double 
      *newLower = getColLower (),
      *newUpper = getColUpper ();

    t_chg_bounds *chgd = new t_chg_bounds [ncols];

    for (int i=0; i<ncols; i++) {
      if (newLower [i] > oldLower [i] + COUENNE_EPS) chgd [i].setLower (t_chg_bounds::CHANGED);
      if (newUpper [i] < oldUpper [i] - COUENNE_EPS) chgd [i].setUpper (t_chg_bounds::CHANGED);
    }

    cutgen_ -> Problem () -> domain () -> push (ncols, NULL, newLower, newUpper);

    if (!(cutgen_ -> Problem () -> btCore (chgd))) // infeasible
      nTightened = -1;

    else {

      const double 
	*newerLower = cutgen_ -> Problem () -> Lb (),
	*newerUpper = cutgen_ -> Problem () -> Ub ();

      for (int i=0; i<ncols; i++) {

	if (newerLower [i] > newLower [i] + COUENNE_EPS) {
	  setColLower (i, newerLower [i]);
	  if (newLower [i] < oldLower [i] + COUENNE_EPS) nTightened++; // extra tightening
	}

      	if (newerUpper [i] < newUpper [i] - COUENNE_EPS) {
	  setColUpper (i, newerUpper [i]);
	  if (newUpper [i] > oldUpper [i] - COUENNE_EPS) nTightened++; // extra tightening
	}
      }
    }

    cutgen_ -> Problem () -> domain () -> pop ();

    delete [] chgd;
  }

  delete [] oldLower;
  delete [] oldUpper;

  return nTightened;
}
