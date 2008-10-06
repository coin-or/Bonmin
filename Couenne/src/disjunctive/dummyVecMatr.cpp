/*
 * Name:    dummyVecMatr.cpp
 * Author:  Pietro Belotti
 * Purpose: fill in empty or single valued vectors and matrices
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CouennePrecisions.hpp"


// take columns of matrix and add each to arrays for matrix under construction
void addSubMatr (int *start, int *len, int *ind, double *el,
		 CoinPackedMatrix &A, 
		 CoinPackedVector &v, 
		 int &cur,
		 int &ncols,
		 int dispM, 
		 int dispVec, 
		 int finalrow) {

  const int 
    *aLe  = A.getVectorLengths (),
    *aIn  = A.getIndices  (),
    *vIn  = v.getIndices  (),
     aCol = A.getMajorDim ();

  int vNum = v.getNumElements ();

  const double 
    *aEl = A.getElements (),
    *vEl = v.getElements ();

  // add each column
  for (int i=0; i<aCol; i++, len++) {

    *start++ = cur;
    *len     = *aLe++;

    // matrix entries
    for (int j = 0; j < *len; j++) {
      *ind++ = dispM + *aIn++;
      *el++  = *aEl++;
    }

    cur += *len;

    // check if there is a corresponding rhs
    if (vNum && (*vIn == i)) {

      ++*len;
      cur++;
      *ind++ = dispVec;
      *el++  = *vEl;

      vIn++;
      vEl++;
      --vNum;
    }

    // normalization entry
    ++*len;
    cur++;
    *ind++ = finalrow;
    *el++  = 1.;

    ++ncols;
  }

  *start = cur;
}


// debug functions
void printMatrix (int nrows, int ncols, int nel, 
		  const int *start, const int *len, 
		  const int *ind,   const double *el) {

  printf ("------------------- %d rows, %d columns, %d nz\n", nrows, ncols, nel);

  for (int i=0, cur = 0; i<nrows; i++) {

    printf ("%2d [%2d -> %2d] (%2d): ", i, start [i], start [i+1] - 1, len [i]);

    for (int j=0; j < len [i]; j++)
      printf ("%d ", ind [start [i] + j]);

    printf (" | --- | ");

    for (int j=0, indice = 0; j < len [i] && j < 1000; j++) {
      while (indice < ind [cur]) {indice++; printf (". ");}
      indice++;
      printf ("%2g ", el [cur++]);
    }

    printf ("\n");
  }
  printf ("-#-\n");
}

void printMatrix (const CoinPackedMatrix *A) {

  int 
    nrows = A -> getMajorDim    (),
    ncols = A -> getMinorDim    (),
    nel   = A -> getNumElements ();

  const double
    *el  = A -> getElements  ();

  const int
    *len   = A -> getVectorLengths (),
    *start = A -> getVectorStarts  (),
    *ind   = A -> getIndices       ();

  printMatrix (nrows, ncols, nel, start, len, ind, el);
}

void printLPMatrix (const OsiSolverInterface &si) {

  // the coefficient matrix
  const CoinPackedMatrix *A = si.getMatrixByCol ();

  printMatrix (A);
}
