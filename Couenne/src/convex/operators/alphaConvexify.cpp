/*
 * Name:    alphaConvexify.cpp
 * Author:  Stefan Vigerske
 * Purpose: create alpha-convexification of a quadratic expression
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>

#include <CoinHelperFunctions.hpp>

#include <OsiSolverInterface.hpp>
#include <IpLapack.hpp>

//#define DEBUG

// fill in one of the two dCoeff vectors
void fill_dCoeff (CouNumber * &, CouNumber, CouNumber *, int);


/** Computes alpha coefficients for an alpha under- and overestimator
 *  of the quadratic term.
 *
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

  if (nqterms_ == 0) {
    nDiag_ = 0;
    return;
  }

  // inverse of dIndex_ mapping: for each variable tell me the
  // index that it will have in dIndex_, or -1 if not there

  int* indexmap = new int [si.getNumCols ()];
  CoinFillN (indexmap, si.getNumCols (), -1);

  if (dIndex_ == NULL)
    make_dIndex (si.getNumCols (), indexmap);
  else
    // build indexmap as inverse of dIndex_
    for (int i=0; i<nDiag_; ++i)
      indexmap [dIndex_ [i]] = i;

  // box diameter
  double *diam = new double [nDiag_];

  const double
    *lower = si.getColLower (),
    *upper = si.getColUpper ();

  for (int i=0; i<nDiag_; ++i) {

    int di = dIndex_ [i];
    diam [i] = upper [di] - lower [di];
  }

  // lower triangular of quadratic term matrix, scaled by box diameter

  double *matrix = new double [nDiag_ * nDiag_];

  CoinFillN (matrix, nDiag_ * nDiag_, 0.);

  for (int i=0; i < nqterms_; ++i) {

    int row = indexmap [qindexI_ [i]];
    int col = indexmap [qindexJ_ [i]];

    // compute value of matrix entry = q_ij * (u_i-l_i) * (u_j-l_j)
    // I (Stefan) do not understand the Lapack docu; it says it needs
    // only the lower triangular but it seem to need both parts to
    // work correct
    double cell = qcoeff_ [i] * diam [row] * diam [col];

    matrix          [col * nDiag_ + row] = cell;
    if (row != col) 
      matrix        [row * nDiag_ + col] = cell;
    //    printf("row %d, col %d: %f\n", row, col, matrix[col*nDiag_+row]);
  }

  delete [] indexmap;

  // compute minimum and maximum eigenvalue of matrix
  double* eigval = new double [nDiag_];
  int info;

  Ipopt::IpLapackDsyev (false,  // do not compute eigenvector
			nDiag_, // dimension
			matrix, // matrix
			nDiag_, // "leading dimension" (number of columns, I think)
			eigval, // output vector to store eigenvalues
			info);  // output status variable
  delete [] matrix;

  if (info != 0) {
    printf ("exprQuad::alphaConvexify: problem computing eigenvalue, info=%d\n", info);
    exit (-1);
    //TODO error handling
  }

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

  delete [] diam;
  delete [] eigval;
}


// fill diagonal vector using eigenvalue passed as parameter 
void fill_dCoeff (CouNumber * &dCoeff, CouNumber eigval, CouNumber *diam, int n) {

  if (dCoeff == NULL)
    dCoeff = new CouNumber [n];

  for (int i=n; i--;) {
    CouNumber di = diam [i];
    dCoeff [i] = (fabs (di) < COUENNE_EPS) ? 
      0. : 
      eigval / (di * di);
  }
}
