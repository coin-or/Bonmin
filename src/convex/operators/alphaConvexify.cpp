/*
 * Name:    alphaConvexify.cpp
 * Author:  Stefan Vigerske
 * Purpose: create alpha-convexification of a quadratic expression
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */


#include "CoinHelperFunctions.hpp"
#include "OsiSolverInterface.hpp"
#include "IpLapack.hpp"

#include "exprQuad.hpp"
#include "CouenneProblem.hpp"

//#define DEBUG

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
 *
 * Return true if a convexification is there, false otherwise.
 */

bool exprQuad::alphaConvexify (const CouenneProblem *p,
			       const OsiSolverInterface &si) {

  if (matrix_.size () == 0)
    return false;

  // inverse of dIndex_ mapping: for each variable tell me the index
  // that it will have in dIndex_, or -1 if not there

  int k=0,
     nDiag    = bounds_.size (),
    *indexmap = new int [si.getNumCols ()],
    *indices  = new int [nDiag];

  CoinFillN (indexmap, si.getNumCols (), -1);

  // box diameter
  double *diam = new double [nDiag];

  bool changed_bounds = false;

  for (std::map <exprVar *, std::pair <CouNumber, CouNumber> >::iterator i = bounds_.begin ();
       i != bounds_.end (); ++i, k++) {

#ifdef DEBUG
  printf ("b%04d. [%20g, %20g]\n", i->first->Index(), i->second.first, i->second.second);
#endif

    int index = i -> first -> Index ();
    indexmap [index] = k;
    indices [k] = index;

    CouNumber
      lb = i -> first -> lb (),
      ub = i -> first -> ub ();

    // if one variable unbounded, bail out
    if ((lb < -COUENNE_INFINITY) ||
	(ub >  COUENNE_INFINITY)) {

      delete [] diam;
      delete [] indexmap;
      delete [] indices;

      return false;
    }

    // if no variable has changed bounds, no need to convexify
    if (fabs (lb - i->second.first)  > COUENNE_EPS) {i -> second.first  = lb; changed_bounds = true;}
    if (fabs (ub - i->second.second) > COUENNE_EPS) {i -> second.second = ub; changed_bounds = true;}

    diam [k] = ub - lb;
#ifdef DEBUG
    printf ("diam %4d - %4d = %g - %g = %g\n", index, k, ub, lb, diam [k]);
#endif
  }

  if (!changed_bounds) {

    delete [] diam;
    delete [] indexmap;
    delete [] indices;

    return true;
  }

  // lower triangular of quadratic term matrix, scaled by box diameter

  double *matrix = new double [nDiag * nDiag];

  CoinFillN (matrix, nDiag * nDiag, 0.);

  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    int 
      xind = row -> first -> Index (),
      irow = indexmap [xind];

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      int 
	yind = col -> first -> Index (),
	icol = indexmap [yind];

      double cell = col -> second * diam [irow] * diam [icol];

      matrix          [icol * nDiag + irow] = cell;
      if (irow != icol) 
	matrix        [irow * nDiag + icol] = cell;
    }
  }

  // compute value of matrix entry = q_ij * (u_i-l_i) * (u_j-l_j)
  // I (Stefan) do not understand the Lapack docu; it says it needs
  // only the lower triangular but it seem to need both parts to
  // work correct

  delete [] indexmap;

  // compute minimum and maximum eigenvalue of matrix
  double* eigval = new double [nDiag];
  int info;

#ifdef DEBUG
  printf ("nDiag = %d\n", nDiag);
  for (int i=0; i<nDiag; i++) {
    for (int j=0; j<nDiag; j++)
      printf ("%6.2f ", matrix [i*nDiag + j]);
    printf ("\n");
  }
#endif 

  Ipopt::IpLapackDsyev (true,   // compute eigenvector
			nDiag,  // dimension
			matrix, // matrix
			nDiag,  // "leading dimension" (number of columns, I think)
			eigval, // output vector to store eigenvalues
			info);  // output status variable

  if (info != 0) {
    printf ("exprQuad::alphaConvexify, warning: problem computing eigenvalue, info=%d\n", info);
    return false;
    //TODO error handling
  }

  // clean eigenvector structure
  eigen_.erase (eigen_.begin (), eigen_.end ());

  for (int i=0; i<nDiag; i++) {

    std::pair <CouNumber, std::vector <std::pair <exprVar *, CouNumber> > > eigenCoord;

    eigenCoord. first = eigval [i];

    for (int j=0; j<nDiag; j++) {

      CouNumber elem = matrix [i * nDiag + j];

      if (fabs (elem) > COUENNE_EPS) 
	eigenCoord. second. push_back (std::pair <exprVar *, CouNumber> 
				       (p -> Var (indices [j]), elem));
    }

    eigen_.push_back (eigenCoord);
  }

#ifdef DEBUG
  for (std::vector <std::pair <CouNumber, 
	 std::vector <std::pair <exprVar *, CouNumber> > > >::iterator i = eigen_.begin ();
       i != eigen_.end (); ++i) {
    printf (" [%g] -- ", i -> first);
    for (std::vector <std::pair <exprVar *, CouNumber> >::iterator j = i -> second. begin();
	 j != i -> second. end (); ++j)
      printf ("(%d,%g) ", j -> first -> Index (), j -> second);
    printf ("\n");
  }
#endif

  delete [] indices;
  delete [] matrix;
  delete [] diam;
  delete [] eigval;

  return true;
}
