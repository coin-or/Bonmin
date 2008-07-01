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

//#define DEBUG

/// generate OsiColCuts for improved (implied and propagated) bounds
void CouenneCutGenerator::genColCuts (const OsiSolverInterface &si,
				      OsiCuts &cs,
				      int nchanged,
				      int *changed) const {

#ifdef DEBUG
  int nrc = cs.sizeRowCuts ();// Must go
#endif

  int  ncols  = problem_ -> nVars (),
      *indLow = new int [ncols], // indices for OsiColCut
      *indUpp = new int [ncols], //
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

#ifdef DEBUG
  for (int i=0; i < problem_ -> nVars (); i++)
    if ((newLow [i] > oldLow [i] + COUENNE_EPS) ||
	(newUpp [i] < oldUpp [i] - COUENNE_EPS))
      printf ("x%-3d. [%-10g , %10g] ---> [%-10g , %10g]\n",
      i, oldLow [i], oldUpp [i], newLow [i], newUpp [i]);*/
#endif

  // check all changed bounds
  for (int i = 0; i < nchanged; i++) {

    int index = changed [i];

    // fails with spectra2 with (abt=2,obbt=0) for variable x70
    //assert (problem_ -> Var (index) -> Multiplicity () > 0);

    if ((index == ind_obj) || 
	(problem_ -> Var (index) -> Multiplicity () <= 0))
      continue;

    if (newLow [index] > newUpp [index])
      problem_ -> Lb (index) = problem_ -> Ub (index);

    CouNumber bd;

    if ((((bd = newLow [index]) > oldLow [index] + COUENNE_EPS) || firstcall_) // better lb?
	&& (bd > -COUENNE_INFINITY / 10)) {                                    // finite?

      //printf ("chging low %d %g -> %g\n", index, oldLow [index], newLow [index]);
      if (problem_ -> Var (index) -> isInteger ()) 
	bd = ceil (bd);
      indLow [nLow]   = index;
      bndLow [nLow++] = bd;
    }

    if ((((bd = newUpp [index]) < oldUpp [index] - COUENNE_EPS) || firstcall_) // better ub?
	&& (bd < COUENNE_INFINITY / 10)) {                                     // finite?

      //printf ("chging upp %d %g -> %g\n", index, oldUpp [index], newUpp [index]);
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

#ifdef DEBUG
  printf ("column cuts\n");
  for (int jj = nrc; jj < cs.sizeRowCuts (); jj++) cs.rowCutPtr (jj) -> print ();
#endif

  delete [] bndLow; delete [] indLow;
  delete [] bndUpp; delete [] indUpp;
}
