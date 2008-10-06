/*
 * Name:    OsiLinear2MatrVec.cpp
 * Author:  Pietro Belotti
 * Purpose: turn OsiSolverInterface objects into coefficient matrix and rhs vector
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"


// construct reduced, standard form matrix M from coefficient matrix of si
void CouenneDisjCuts::OsiSI2MatrVec (CoinPackedMatrix &M, 
				     CoinPackedVector &r,
				     OsiSolverInterface &si) const {

  // the coefficient matrix
  const CoinPackedMatrix *A = si.getMatrixByRow ();

  int 
    nrows = A -> getMajorDim    (),
    ncols = A -> getMinorDim    (),
    nel   = A -> getNumElements ();

  const char *sen = si.getRowSense ();

  const double
    *rac = si.getRowActivity (),
    *x   = si.getColSolution (),
    *rlb = si.getRowLower (),
    *rub = si.getRowUpper (),
    *clb = si.getColLower (),
    *cub = si.getColUpper (),
    *el  = A -> getElements  ();

  const int
    *len   = A -> getVectorLengths (),
    *start = A -> getVectorStarts  (),
    *ind   = A -> getIndices       ();

  /////////////////// count rows that will have to be duplicated
  int 
    ndupEl   = 0, // count nonzero elements
    ndupRows = 0; // count rows

  for (int i=0; i<nrows; i++, sen++, len++) 
    if ((*sen == 'E') || 
	(*sen == 'R')) {
      ndupEl += *len;
      ndupRows++;
    }

  sen -= nrows;
  len -= nrows;

  double *mEl = new double [nel + ndupEl + 2*ncols];

  r.reserve (nrows + ndupRows + 2*ncols + 1);

  int 
    mRows = 0,
    curA  = 0, 
    curM  = 0, 
    *mIn  = new int [nel + ndupEl + 2*ncols],
    *mSt  = new int [nrows + ndupRows + 2*ncols + 1],
    *mLe  = new int [nrows + ndupRows + 2*ncols];


  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("matrix A (%d %d) %d elements:\n", nrows, ncols, nel);

    printf ("start: ");   for (int i=0; i <= nrows; i++) printf ("%d ", start [i]);
    printf ("\nlen:   "); for (int i=0; i <  nrows; i++) printf ("%d ", len   [i]);

    printf ("\nElements:\n"); 
    for (int i=0, j=0; i<nrows; i++) {
      for (int k=0; k<start [i+1] - start [i]; k++, j++) 
	printf ("(%d %g) ", ind [j], el [j]);
      printf (" in [%g,%g]\n", rlb [i], rub [i]);
    }
  }

  // for each constraint
  //   if activerows on and no activity, continue
  //   if 'E' or 'L' or 'R', copy    entries in matrix and    upper bound in R
  //   if 'E' or 'G' or 'R', copy -1*entries in matrix and -1*lower bound in R
  //
  // for each variable
  //   if activecols on and within bounds, continue
  //   if upper bounded, copy  1 in matrix and    upper bound in R
  //   if lower bounded, copy -1 in matrix and -1*lower bound in R

  // Constraints //////////////////////////////////
  for (int i=0; i<nrows; i++, rac++, rlb++, rub++, sen++, start++, curA += *len++) {

    if (activeRows_ &&
	(*rac < *rub - COUENNE_EPS) &&
	(*rac > *rlb + COUENNE_EPS))
      continue;

    if (*sen != 'G') {
      *mSt++ = curM;
      *mLe++ = *len;
      CoinCopyN (ind + curA, *len, mIn + curM);
      CoinCopyN (el  + curA, *len, mEl + curM);
      curM += *len;
      if (*rub != 0.) r.insert (mRows, *rub);
      mRows++;
    }

    if (*sen != 'L') {
      *mSt++ = curM;
      *mLe++ = *len;
      CoinCopyN (ind + curA, *len, mIn + curM);
      CoinInvN  (el  + curA, *len, mEl + curM); // invert contents
      curM += *len;
      if (*rlb != 0.) r.insert (mRows, -*rlb); // invert bound
      mRows++;
    }
  }


  // Variables ////////////////////////////////////
  for (int i=0; i<ncols; i++, x++, clb++, cub++) {

    if (activeCols_ && 
	(*x < *cub - COUENNE_EPS) &&
	(*x > *clb + COUENNE_EPS) ||
	(couenneCG_ -> Problem () -> Var (i) -> Multiplicity () <= 0))
      continue;

    if (*clb > -COUENNE_INFINITY) { // has lower bound -- invert!
      *mSt++ = curM;
      *mLe++ = 1;
      mIn [curM]   = i;
      mEl [curM++] = -1.;
      if (*clb != 0.) r.insert (mRows, -*clb);  // invert bound so x >= b becomes -x <= -b
      mRows++;
    }

    if (*cub <  COUENNE_INFINITY) { // has upper bound
      *mSt++ = curM;
      *mLe++ = 1;
      mIn [curM] = i;
      mEl [curM++] = 1.;
      if (*cub != 0.) r.insert (mRows, *cub);
      mRows++;
    }
  }

  *mSt = curM;

  mSt -= mRows;
  mLe -= mRows;

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("matrix (%d %d) %d elements:\n", mRows, ncols, curM);

    printf ("start: ");   for (int i=0; i<=mRows; i++) printf ("%d ", mSt [i]);
    printf ("\nlen:   "); for (int i=0; i<mRows;  i++) printf ("%d ", mLe [i]);

    printf ("\nElements:\n"); 
    for (int i=0, j=0; i<mRows; i++) {
      for (int k=0; k<mSt[i+1] - mSt[i]; k++, j++) 
	printf ("(%d %g) ", mIn [j], mEl [j]);
      printf ("\n");
    }

    printf ("vector:\n"); 
    for (int i=0; i<r.getNumElements (); i++)
      printf ("(%d %g) ", r.getIndices () [i], r.getElements () [i]);
    printf ("\n");
  }

  M.assignMatrix (true,     // column ordered
		  ncols,    // minor dimension
		  mRows,    // major dimension
		  curM,     // number of elements
		  mEl,      // elements
		  mIn,      // indices
		  mSt,      // starting positions
		  mLe);     // length
}
