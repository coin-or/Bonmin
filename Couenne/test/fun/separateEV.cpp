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

  double *A = (double *) malloc (np*np * sizeof (double));

  {
    register int i,j;

    A [0] = 1;

    for (j=1; j<np; j++)
      A [np*j] = sol [j-1];

    for (i=1; i<np; i++)
      for (j=i; j<np; j++)
	A [np*j+i] = sol [indexQ (i-1,j-1,n)];
  }

  int m;
  double *w, *z;

  dsyevx_wrapper (    n, A, m, w, z);
  eigenPlay      (cs, n,    m, w, z);

  free (A); 

  ncalls++;
}


///
void SdpCutGen::dsyevx_wrapper (int n, double *A, int &m, double *w, double *z) const {

  w = (double *) malloc (n   * sizeof (double));
  z = (double *) malloc (n*n * sizeof (double));

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
    //    m      = n,   // output number of eigenvalues found

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
  
  if (info!=0) 
    printf ("::: dsyevx returned status %d\n", info);

  if ((info==0) && (m > 0)) { // there is at least one eigenvector

    printf ("%d negative eigenvalues found: ", m);

    for (int i=0; i<m; i++) 
      printf ("%.3f ", w [i]);
    printf ("\n");

    //    eigenPlay (n, m, w, z);
  }

  free (work);
  free (ifail);
  free (iwork);
}
