/*
 * Name:    singleDisjunctions.cpp
 * Author:  Pietro Belotti
 * Purpose: simpler methods on single disjunctions
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "CouenneDisjCuts.hpp"


/// create single osicolcut disjunction
OsiCuts *CouenneDisjCuts::getSingleDisjunction (OsiSolverInterface &si) const { 

  int 
    ncols = si.getNumCols (), nNL = 0, nNU = 0,
    *indNL = new int [ncols],      
    *indNU = new int [ncols];

  double 
    *oldL  = problem_ -> Lb (),   *valNL = new double [ncols],
    *oldU  = problem_ -> Ub (),   *valNU = new double [ncols];

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

  //cut. print ();

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

  printf ("into getBoxUnion\n");

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
    lowerLeft. append (right -> colCutPtr (i) -> lbs ());
    upperLeft. append (right -> colCutPtr (i) -> ubs ());
  }

  // sort all indexed vectors
  lowerLeft.  sortIncrIndex ();  upperLeft.  sortIncrIndex ();
  lowerRight. sortIncrIndex ();  upperRight. sortIncrIndex ();

  // now compare bounds one by one -- complexity linear in size of vectors

  double newBound;

  {   // lower bounds /////////////////////////////////////////////////

    const register int 
      *lLi = lowerLeft.  getIndices (),
      *lRi = lowerRight. getIndices ();

    const double
      *lLe = lowerLeft.  getElements (),
      *lRe = lowerRight. getElements ();
      //*oldLower = si.getColLower ();

    int
      lLn = lowerLeft.  getNumElements (),
      lRn = lowerRight. getNumElements ();

    for (;;) {

      for (;;) {
	register int diff = *lLi - *lRi;
	if      (diff < 0) {if (!--lLn) break; lLi++; lLe++;}
	else if (diff > 0) {if (!--lRn) break; lRi++; lRe++;}
	else break;
      }

      if (!lLn || !lRn) break;                                    // here, !lLn, !lRn, or
      lower. insert (*lLi, newBound = *lLe < *lRe ? *lLe : *lRe); // *lLi==*lRi: add min(*lLe, *lRe)) 

      printf ("common bound: x_%d >= %g\n", *lLi, newBound);

      // check if tighten problem
      //if (newBound > oldLower [*lLi] + COUENNE_EPS) retval = COUENNE_TIGHTENED;
    }
  }

  {  // upper bounds /////////////////////////////////////////////////

    const register int 
      *uLi = upperLeft.  getIndices (),
      *uRi = upperRight. getIndices ();

    const double
      *uLe = upperLeft.  getElements (),
      *uRe = upperRight. getElements ();
      //*oldUpper = si.getColUpper ();

    int
      uLn = upperLeft.  getNumElements (),
      uRn = upperRight. getNumElements ();

    for (;;) {

      for (;;) {
	register int diff = *uLi - *uRi;
	if      (diff < 0) {if (!--uLn) break; uLi++; uLe++;}
	else if (diff > 0) {if (!--uRn) break; uRi++; uRe++;}
	else break;
      }

      if (!uLn || !uRn) break;                                    // here, !uLn, !uRn, or
      upper. insert (*uLi, newBound = *uLe < *uRe ? *uLe : *uRe); // *uLi==*uRi: add min(*uLe, *uRe)) 

      printf ("common bound: x_%d <= %g\n", *uLi, newBound);

      // check if tighten problem
      //if (newBound < oldUpper [*uLi] - COUENNE_EPS) retval = COUENNE_TIGHTENED;
    }
  }

  printf ("============\n");

  return retval;
}
