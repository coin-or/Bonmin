#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>


extern "C" {

  /* Lapack routine to get an LU decomposition of matrix A */
  void dgetrf_ (int    *,  /* M: number of rows    of A    */
		int    *,  /* N:           columns         */
		double *,  /* A: the matrix in column ord. */
		int    *,  /* LDA: leading dim. of A       */
		int    *,  /* IPIV: permutation data       */
		int    *); /* INFO: result                 */

  /* Compute inverse of a triangular matrix                 */
  void dtptri_ (char   *,  /* UPLO: upper or lower?         */
		char   *,  /* DIAG: only 1 on diagonal?     */
		int    *,  /* N:    order                   */
		double *,  /* AP:   A (packed)              */
		int    *); /* INFO: output status           */
}


#define MAX_OUT_ITER 20
#define MAX_EIGENVALUE -1e-6

// sdpcut separator
void SdpCutGen::separateLU (const OsiSolverInterface &si, 
			    OsiCuts &cs) const {
  /*
  int n  = *(const_cast <int *> (&n_)), 
      np = n + 1;

  const double *sol = si.getColSolution ();
  */
  /* 
   *  get LU decomposition of 
   *
   *          / 1  x^T \
   *  X'  =   \ x  X   /
   *
   *  obtaining X_0 = PLU
   *
   *  as U is upper triangular and contains positive and negative
   *  pivots on the diagonal, the inequality e_j^T U e_j < 0 is
   *  violated for all j such that U_{jj} < 0.
   *  
   *  Therefore, as L^{-1} P^T X_0 = U we have
   *
   *  e_j^T L^{-1} P^T X_0 e_j >= 0 is a valid inequality
   */
  /*
  int 
    nrows = np,                                 // # rows    of A
    ncols = np,                                 //   columns
    *ipiv = (int *) malloc (np * sizeof (int)), // permutation data
    info,                                       // output status
    lda   = np;                                 // leading dimension of A

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

  nrows = 3;
  ncols = 3;
  A [0] =   0; A [3] = 9; A [6] =  0;
  A [1] =   0; A [4] = 0; A [7] =  4;
  A [2] =  -5; A [5] = 0; A [8] =  0;
  lda = 3;

  dgetrf_ (&nrows, &ncols, A, &lda, ipiv, &info);
  */
  /*
   *  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
   *          On entry, the upper or lower triangular matrix A, stored
   *          columnwise in a linear array.  The j-th column of A is stored
   *          in the array AP as follows:
   *          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
   *          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.
   *          See below for further details.
   *          On exit, the (triangular) inverse of the original matrix, in
   *          the same packed storage format.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular
   *                matrix is singular and its inverse can not be computed.
   *
   */
  /*
  printf ("[L U] =\n");

  for (int i=0; i<np; i++) {
    for (int j=0; j<=i; j++) 
      printf ("%+6.2f ", A [np*j+i]);
    printf ("                        ");
    for (int j=i+1; j<np; j++) 
      printf ("%+6.2f ", A [np*j+i]);
    printf ("\n");
  }

  for (int i=0; i<np; i++)
    printf ("%d --> %d\n", i, ipiv [i] - 1);

  for (int i=0; i<np; i++) {
    for (int j=0; j<np; j++) 
      printf ("%+4.1f  ", A [3*j+i]);
    printf ("\n");
  }

  int *P  = (int *) malloc (np * sizeof (int));
  int *PT = (int *) malloc (np * sizeof (int));

  for (int i=0, *pp = ipiv; i<np; i++) {
    if (*pp == i);    
  }

  int N = (n+1) * (n+2) / 2; // size of packed lower triangular matrix L to invert

  double *AP = (double *) malloc (N * sizeof (double));

  for (int i=0; i<np; i++)
    for (int j=0; j<=i; j++) 
      AP [i + (j-1)*((2*np - j)/2)] = A [np * j + i];

  char uplo = 'L';
  char diag = 'U';

  dtptri_ (&uplo, &diag, &np, AP, &info);

  if (info)
    printf ("::: SeparateLU: dtptri returned nonzero\n");



  for (int i=0; i<np; i++) {
    for (int j=0; j<np; j++)
      if (A [3*j+i] != A [3*i+j]) printf ("!");
      else if ((i==j) && (A [3*j+i] < 0)) printf ("<");
  }
  printf ("\n");

  if (info==0) { // there is at least one eigenvector

    int N = n*(n+3)/2;

    double *x = (double *) malloc (np * sizeof (double));

    free (x);
  }

  free (ipiv);
  free (A);
  */

}
