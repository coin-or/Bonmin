/*
 * Name:    singleDisjunctions.cpp
 * Author:  Pietro Belotti
 * Purpose: simpler methods on single disjunctions
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneDisjCuts.hpp"


/// create single osicolcut disjunction
OsiCuts *CouenneDisjCuts::getSingleDisjunction (OsiSolverInterface &si) const { 

  int 
    ncols = si.getNumCols (), nNL = 0, nNU = 0,
    *indNL = new int [ncols],      
    *indNU = new int [ncols];

  double 
    *oldL  = couenneCG_ -> Problem () -> Lb (),   *valNL = new double [ncols],
    *oldU  = couenneCG_ -> Problem () -> Ub (),   *valNU = new double [ncols];

  const double
    *newL = si.getColLower (),   
    *newU = si.getColUpper ();

  for (int i=0; i<ncols; i++) {
    if (newL [i] > oldL [i] + COUENNE_EPS) {indNL [nNL] = i; valNL [nNL++] = newL [i];}
    if (newU [i] < oldU [i] - COUENNE_EPS) {indNU [nNU] = i; valNU [nNU++] = newU [i];}
  }	      

  OsiColCut cut;

  cut. setLbs (nNL, indNL, valNL);
  cut. setUbs (nNU, indNU, valNU);

  OsiCuts *cuts = new OsiCuts;

  cuts -> insert (cut);

  delete [] indNL;   delete [] valNL;
  delete [] indNU;   delete [] valNU;

  return cuts;
}


/// check if (column!) cuts compatible with solver interface
int CouenneDisjCuts::checkDisjSide (OsiSolverInterface &si, OsiCuts *cuts) const {

  int retval = COUENNE_FEASIBLE;

  const double 
    *lower = si.getColLower (),
    *upper = si.getColUpper ();

  // check each colcut in cuts (there should be just one)

  for (int i = cuts -> sizeColCuts (); i--;) {

    // lower bounds

    const CoinPackedVector &lbs = cuts -> colCutPtr (i) -> lbs ();
    const int    *lindices = lbs.getIndices ();
    const double *lvalues  = lbs.getElements ();

    for (int j = lbs.getNumElements (); j--;) {
      register double lb  = *lvalues++;
      register int    ind = *lindices++;

      if (lb > upper [ind] + COUENNE_EPS) // fathom node
	return COUENNE_INFEASIBLE;

      if (lb > lower [ind] + COUENNE_EPS) 
	retval = COUENNE_TIGHTENED;
    }

    // upper bounds

    const CoinPackedVector &ubs = cuts -> colCutPtr (i) -> ubs ();
    const int    *uindices = ubs.getIndices ();
    const double *uvalues  = ubs.getElements ();

    for (int j = ubs.getNumElements (); j--;) {
      register double ub  = *uvalues++;
      register int    ind = *uindices++;

      if (ub < lower [ind] - COUENNE_EPS) // fathom node
	return COUENNE_INFEASIBLE;

      if (ub < upper [ind] - COUENNE_EPS) 
	retval = COUENNE_TIGHTENED;
    }
  }

  return retval;
}


/// compute smallest box containing both left and right boxes.
int CouenneDisjCuts::getBoxUnion (OsiSolverInterface &si,
				  OsiCuts *left, OsiCuts *right, 
				  CoinPackedVector &lower, CoinPackedVector &upper) const {

  int retval = COUENNE_FEASIBLE;

  CoinPackedVector 
    lowerLeft,  upperLeft,
    lowerRight, upperRight;

  // merge all left lowers and uppers
  for (int i = left -> sizeColCuts (); i--;) {
    lowerLeft. append (left -> colCutPtr (i) -> lbs ());
    upperLeft. append (left -> colCutPtr (i) -> ubs ());
  }

  // merge all right lowers and uppers
  for (int i = right -> sizeColCuts (); i--;) {
    lowerRight. append (right -> colCutPtr (i) -> lbs ());
    upperRight. append (right -> colCutPtr (i) -> ubs ());
  }

  // sort all indexed vectors
  lowerLeft.  sortIncrIndex ();  upperLeft.  sortIncrIndex ();
  lowerRight. sortIncrIndex ();  upperRight. sortIncrIndex ();

  // now compare bounds one by one -- complexity linear in size of vectors

  mergeBoxes (-1, lowerLeft, lowerRight, lower);
  mergeBoxes (+1, upperLeft, upperRight, upper);

  return retval;
}


/// utility to merge vectors into one
void CouenneDisjCuts::mergeBoxes (int dir, // direction (negative for "<", positive for ">")
				  CoinPackedVector &left, 
				  CoinPackedVector &right, 
				  CoinPackedVector merged) const { 
  int
    Ln = left.  getNumElements (),
    Rn = right. getNumElements ();

  if (!Ln || !Rn) 
    return;

  const int 
    *Li = left.  getIndices (),
    *Ri = right. getIndices ();

  const double
    *Le = left.  getElements (),
    *Re = right. getElements ();

  for (;;) {

    for (;;) {

      register int diff = *Li - *Ri;

      if      (diff < 0) {if (!--Ln) break; Li++; Le++;}
      else if (diff > 0) {if (!--Rn) break; Ri++; Re++;}
      else break;
    }

    if (!Ln || !Rn) break;                                  // !Ln, !Rn (==> exit), or Li==*Ri: 
    if (dir < 0) merged. insert (*Li, *Le<*Re ? *Le : *Re); // add min(*Le, *Re)) 
    else         merged. insert (*Li, *Le>*Re ? *Le : *Re); // add max(*Le, *Re)) 

    Li++; Ri++;
    Le++; Re++;

    if (!--Ln || !--Rn) break;
  }
}
