/*
 * Name:    genColCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate Column Cuts for improved bounds
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

/// generate OsiColCuts for improved (implied and propagated) bounds
void CouenneCutGenerator::genColCuts (const OsiSolverInterface &si,
				      OsiCuts &cs,
				      int nchanged, 
				      int *changed) const {

  int ncols = problem_ -> nVars () + problem_ -> nAuxs (),
    *indLow = new int [ncols], // indices for OsiColCut
    *indUpp = new int [ncols],
    nLow, nUpp = nLow = 0;

  // values fo OsiColCut
  CouNumber *bndLow = new CouNumber [ncols],
    *bndUpp = new CouNumber [ncols];

  const CouNumber 
    *oldLow = si.getColLower (), // old bounds
    *oldUpp = si.getColUpper (),
    *newLow = problem_ -> Lb (), // changed bounds
    *newUpp = problem_ -> Ub ();

  // check all changed bounds
  for (int i=0; i<nchanged; i++) {

    register int index = changed [i];

    CouNumber bd;

    if (firstcall_ || ((bd = newLow [index]) > oldLow [index] + COUENNE_EPS)) { // lower
      indLow [nLow]   = index;
      bndLow [nLow++] = bd;
    }

    if (firstcall_ || ((bd = newUpp [index]) < oldUpp [index] - COUENNE_EPS)) { // upper
      indUpp [nUpp]   = index;
      bndUpp [nUpp++] = bd;
    }
  }

  // create Column Cut

  if (nUpp || nLow) {

    OsiColCut *cut = new OsiColCut;

    if (cut) {
      cut -> setLbs (nLow, indLow, bndLow);
      cut -> setUbs (nUpp, indUpp, bndUpp);

      cs.insert (cut);
      delete cut;
    }
  }

  delete [] bndLow; delete [] indLow;
  delete [] bndUpp; delete [] indUpp;
}
