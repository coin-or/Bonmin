/*
 * Name:    alphaConvexify.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: create alpha-convexification of a quadratic expression
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>

#include <OsiSolverInterface.hpp>
#include <IpLapack.hpp>


// fill in one of the two dCoeff vectors
void fill_dCoeff (CouNumber * &, CouNumber, CouNumber *, int);

/** Computes alpha coefficients for an alpha under- and overestimator of the quadratic term.
 * For the underestimator, dCoeffLo_ is computed such that
 *
 * x^TQx + sum_i dCoeffLo_i (x_i - lb x_i)(ub x_i - x-i) is convex and
 * underestimating (alpha_i is negative),
 *
 * Regarding the overestimator, dCoeffUp_ are computed such that
 *
 * x^TQx + sum_i dCoeffUp_i (x_i - lb x_i)(ub x_i - x-i) is concave
 * and overestimating (alpha_i is positive).
 *
 * If the method hasn't been called, dIndex_ will be NULL.
 * If Q is positive semidefinite, then dCoeffLo_ will be NULL.
 * If Q is negative semidefinite, then dCoeffUp_ will be NULL.
*/

void exprQuad::alphaConvexify (const OsiSolverInterface &si) {

  if (getnQTerms () == 0) {
    nDiag_ = 0;
    return;
  }

  // inverse of dIndex_ mapping, for each variable tell me the
  // index that it will have in dIndex_, or -1 if not there

  int* indexmap = new int [si.getNumCols ()];

  for (int i=0; i<si.getNumCols(); ++i)
    indexmap [i] = -1;

  if (dIndex_ == NULL) { 

    // first time called... check which variables are there, and where
    // we will put in the dIndex_ array

    int *qindexI = getQIndexI (),
        *qindexJ = getQIndexJ ();

    nDiag_=0;

    // fill indexmap
    for (int i=0; i<getnQTerms(); ++i) {

      int qi = qindexI [i], 
	  qj = qindexJ [i];

      if                (indexmap [qi] == -1)  indexmap [qi] = nDiag_++;
      if ((qi != qj) && (indexmap [qj] == -1)) indexmap [qj] = nDiag_++;
    }

    dIndex_ = new int [nDiag_];

    for (int i=0; i<si.getNumCols(); ++i) {
      if (indexmap[i] > -1) {
	dIndex_ [indexmap [i]] = i;
      }
    }
  } else {
    // build indexmap as inverse of dIndex_
    for (int i=0; i<nDiag_; ++i)
      indexmap [dIndex_ [i]] = i;
  }

  // box diameter
  double* diam = new double [nDiag_];
  for (int i=0; i<nDiag_; ++i)
    diam [i] = 
      si.getColUpper () [dIndex_ [i]] -
      si.getColLower () [dIndex_ [i]];

  // lower triangular of quadratic term matrix, scaled by box diameter
  double* matrix = new double [nDiag_ * nDiag_];

  for (int i=0; i < nDiag_ * nDiag_; ++i)
    matrix [i] = 0.;

  for (int i=0; i<getnQTerms(); ++i) {

    int row = indexmap [getQIndexI () [i]];
    int col = indexmap [getQIndexJ () [i]];

    // compute value of matrix entry = q_ij * (u_i-l_i) * (u_j-l_j)
    // I (Stefan) do not understand the Lapack docu; it says it needs only the lower triangular
    // but it seem to need both parts to work correct
                    matrix [col * nDiag_ + row] = getQCoeffs () [i] * diam [row] * diam [col];
    if (row != col) matrix [row * nDiag_ + col] = getQCoeffs () [i] * diam [row] * diam [col];
    //		printf("row %d, col %d: %f\n", row, col, matrix[col*nDiag_+row]);
  }

  // compute minimum and maximum eigenvalue of matrix
  // ok, currently computes all eigenvalues
  double* eigval = new double [nDiag_];
  int info;

  Ipopt::IpLapackDsyev (false, nDiag_, matrix, nDiag_, eigval, info);

  if (info != 0) {
    printf ("exprQuad::alphaConvexify: problem computing eigenvalue, info=%d\n", info);
    return;
    //TODO error handling
  }
  //	printf("min. eigvalue: %f\n", eigval[0]);

  // if min. eigenvalue negative, setup dCoeffLo_
  if (eigval [0] < 0) 
    fill_dCoeff (dCoeffLo_, eigval [0], diam, nDiag_);
  else  // quadratic term is convex, no convexification needed
    if (dCoeffLo_) {
      delete dCoeffLo_;
      dCoeffLo_ = NULL;
    }

  // if max. eigenvalue is positive, setup dCoeffUp_
  if (eigval [nDiag_ - 1] > 0)
    fill_dCoeff (dCoeffUp_, eigval [nDiag_ - 1], diam, nDiag_);
  else // quadratic term is concave, no "concavification" needed
    if (dCoeffUp_) {
      delete dCoeffUp_;
      dCoeffUp_ = NULL;
    }

  delete[] matrix;
  delete[] diam;
  delete[] eigval;
  delete[] indexmap;
}


// fill diagonal vector using eigenvalue passed as parameter 
void fill_dCoeff (CouNumber * &dCoeff, CouNumber eigval, CouNumber *diam, int n) {

  if (dCoeff == NULL)
    dCoeff = new CouNumber [n];
  for (int i=0; i<n; ++i) {
    if (diam [i] == 0.) dCoeff [i] = 0.;
    else                dCoeff [i] = eigval / (diam [i] * diam [i]);
  }
}
