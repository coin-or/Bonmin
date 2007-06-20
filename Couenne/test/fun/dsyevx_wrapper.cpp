/*
 * Name:    dsyevx_wrapper.cpp
 * Author:  Pietro Belotti
 * Purpose: wrapper for Lapack's Fortran routine to compute eigen*
 *
 * (C) Pietro Belotti. This file is publised under the Common Public License.
 */

#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>


extern "C" {

  /* Lapack routine to compute orthonormal eigenvalues/eigenvectors (in Fortran) */
  void dsyevx_ (char   *,
		char   *,
		char   *,
		int    *,
		double *,
		int    *,
		double *,
		double *,
		int    *,
		int    *,
		double *,
		int    *,
		double *,
		double *,
		int    *,
		double *,
		int    *,
		int    *,
		int    *,
		int    *);
}


#define MAX_EIGENVALUE 1e-1
//#define MAX_EIGENVALUE COIN_DBL_MAX

///
void dsyevx_wrapper (int n, double *A, int &m, double * &w, double * &z) {

  w = (double *) malloc (n   * sizeof (double));
  z = (double *) malloc (n*n * sizeof (double));
  m = n;

  static int lwork = -1;

  if (lwork < 0) 
    lwork = 8*n;

  char 
    jobz  = 'V',  // compute both eigenvalues and eigenvectors
    range = 'V',  // range for selection is on values of eigenvalues
    uplo  = 'U';  // upper triangular matrix is given

  int
    il     = 1,   // first  eigenvalue to be returned (not used)
    iu     = n,   // second                           (not used)
    info,         // output status
    lda    = n,   // leading dimension of A
    ldz    = n,   // leading dimension of z

    *ifail = (int *) malloc (n   * sizeof (int)),
    *iwork = (int *) malloc (5*n * sizeof (int));
  
  double
    abstol = 1e-7,                      // absolute tolerance
    vl     = -COIN_DBL_MAX,             // minimum eigenvalue wanted
    vu     = MAX_EIGENVALUE,            // maximum
    *work  = (double *) malloc (lwork * sizeof (double));

  dsyevx_ (&jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu,
	   &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

  lwork = (int) (*work);

  // create corresponding cut ////////////////////////////////////
  
  if (info) 
    printf ("::: dsyevx returned status %d\n", info);

  if ((info == 0) && (m > 0)) { // there is at least one eigenvector

    z = (double *) realloc (z, m*n * sizeof (double));
    w = (double *) realloc (w, m   * sizeof (double));

    printf ("%d negative eigenvalues found: ", m);

    for (int i=0; i<m; i++) 
      printf ("%.3f ", w [i]);
    printf ("\n");
  }

  free (work);
  free (ifail);
  free (iwork);
}
