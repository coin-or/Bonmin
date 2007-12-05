/*
 * Name:    genColCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate Column Cuts for improved bounds
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

/// generate OsiColCuts for improved (implied and propagated) bounds
void CouenneCutGenerator::genColCuts (const OsiSolverInterface &si,
				      OsiCuts &cs,
				      int nchanged,
				      int *changed) const {

  int  ncols  = problem_ -> nVars (),
      *indLow = new int [ncols], // indices for OsiColCut
      *indUpp = new int [ncols],
       nLow, nUpp = nLow = 0,
       ind_obj = problem_ -> Obj (0) -> Body () -> Index ();

  // values fo OsiColCut
  CouNumber *bndLow = new CouNumber [ncols],
            *bndUpp = new CouNumber [ncols];

  const CouNumber 
    *oldLow = si.getColLower (), // old bounds
    *oldUpp = si.getColUpper (),
    *newLow = problem_ -> Lb (), // changed bounds
    *newUpp = problem_ -> Ub ();

  // check all changed bounds
  for (int i = 0; i < nchanged; i++) {

    register int index = changed [i];

    if (index == ind_obj)
      continue;

    if (newLow [index] > newUpp [index])
      problem_ -> Lb (index) = problem_ -> Ub (index);

    CouNumber bd;

    if (((bd = newLow [index]) > oldLow [index] + COUENNE_EPS) || firstcall_) { // lower

      if (problem_ -> Var (index) -> isInteger ()) 
	bd = ceil (bd);
      indLow [nLow]   = index;
      bndLow [nLow++] = bd;
    }

    if (((bd = newUpp [index]) < oldUpp [index] - COUENNE_EPS) || firstcall_) { // upper

      if (problem_ -> Var (index) -> isInteger ()) 
	bd = floor (bd);
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
