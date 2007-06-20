#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

extern "C" {

  /* Lapack routine to compute eigenvalues/eigenvectors (in Fortran) */
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


#define MAX_OUT_ITER 20
#define MAX_EIGENVALUE -1e-2
//#define MAX_EIGENVALUE COIN_DBL_MAX

// sdpcut separator
void SdpCutGen::separateEV (const OsiSolverInterface &si, 
			    OsiCuts &cs) const {

  static int ncalls = 0;

  int n  = *(const_cast <int *> (&n_)), 
      np = n + 1;

  const double *sol = si.getColSolution ();

  /*
   *  get eigenvectors v for all negative eigenvalues of
   *
   *         / 1  x' \
   *  X' =   \ x  X  /
   *
   *  And use them to add the valid inequality 
   *
   *  v^T X' v >= 0
   */

  static int lwork = -1;

  if (lwork < 0) 
    lwork = 8*np;

  if (ncalls < 5000)  {

    char 
      jobz  = 'V',   // compute both eigenvalues and eigenvectors
      range = 'V',   // range for selection is on values of eigenvalues
      uplo  = 'U';   // upper triangular matrix is given

    int
      il     = 1,    // first  eigenvalue to be returned (not used)
      iu     = np,   // second                           (not used)
      info,          // output status
      lda    = np,   // leading dimension of A
      ldz    = np,   // leading dimension of z
      m      = np,   // output number of eigenvalues found

      *ifail = (int *) malloc (np   * sizeof (int)),
      *iwork = (int *) malloc (5*np * sizeof (int));

    double
      abstol = 1e-7,                      // absolute tolerance
      vl     = -COIN_DBL_MAX,             // minimum eigenvalue wanted
      vu     = MAX_EIGENVALUE,            // maximum
      *A     = (double *) malloc (np*np * sizeof (double)),
      *w     = (double *) malloc (np    * sizeof (double)),
      *work  = (double *) malloc (lwork * sizeof (double)),
      *z     = (double *) malloc (np*np * sizeof (double));

    {
      register int i,j;

      A [0] = 1;

      for (j=1; j<np; j++)
	A [np*j] = sol [j-1];

      for (i=1; i<np; i++)
	for (j=i; j<np; j++)
	  A [np*j+i] = sol [indexQ (i-1,j-1,n)];
    }

    dsyevx_ (&jobz, &range, &uplo, &np, A, &lda, &vl, &vu, &il, &iu,
	     &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

    lwork = (int) (*work);

    // create corresponding cut ////////////////////////////////////

    if (info!=0) 
      printf ("::: dsyevx returned status %d\n", info);

    if ((info==0) && (m > 0)) { // there is at least one eigenvector

      printf ("%d negative eigenvalues found: ", m);

      for (int i=0; i<m; i++) 
	printf ("%.3f ", w [i]);
      printf ("\n");

      double *x = (double *) malloc (np * sizeof (double));

      for (int k=0; k<m; k++) { // eigenvalue negative enough

	/*
	printf ("%3d: %6.3f -- ", k, w [k]);
	for (int j=0; j<np; j++) 
	  printf ("%.3f ", z [k*np + j]);
	printf ("\n");
	*/

	double *zbase = z + k*np;

#define EIG_TOL 0

	for (int j=0; j<np; j++) {
	  x [j] = *zbase++;
	  if (fabs (x [j]) < EIG_TOL) 
	    x [j] = 0;
	}

	genSDPcut (cs, x, x);
      }

      free (x);
    }

    free (A); free (w);
    free (z); free (work);
    free (ifail);
    free (iwork);
  }

  ncalls++;
}
