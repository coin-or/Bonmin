#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

extern "C" {
  /* Fortran routine solving the TRP */
  void dgqt_ (int    *,  /* order of matrix A            */
	      double *,  /* matrix A                     */
	      int    *,  /* leading dimension of A       */
	      double *,  /* vector b                     */
	      double *,  /* rhs of constraint            */
	      double *,  /* rel. tol.                    */
	      double *,  /* abs. tol.                    */
	      int    *,  /* max iterations               */
	      double *,  /* estimate multiplier of con.  */
	      double *,  /* objective function value     */
	      double *,  /* vector x                     */
	      int    *,  /* output status                */
	      double *,  /* vector                       */
	      double *,  /* vector                       */
	      double *); /* vector                       */

  /* Fortran routine to obtain eigen{values,vectors} */
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
#define MAX_EIGENVALUE -1e-6

// sdpcut separator
void SdpCutGen::generateCuts (const OsiSolverInterface &si, 
			      OsiCuts &cs, 
			      const CglTreeInfo info) const {
  static int ncalls = 0;

  int n  = *(const_cast <int *> (&n_)), 
      np = n + 1;

  int N = n*(n+3)/2;

  const double *sol = si.getColSolution ();

  /* TODO: 
   *
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

  if (ncalls < 50)  {

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
      *w     = (double *) malloc (n     * sizeof (double)),
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

    if ((info==0) && (m > 1)) { // there is at least one eigenvector

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

	OsiRowCut *cut = new OsiRowCut;

	double *coeff  = new double [N];
	int    *ind    = new int    [N],
   	        nterms = 0;

	double *zbase = z + k*np;

#define EIG_TOL 0

	for (int j=0; j<np; j++) {
	  x [j] = *z++;
	  if (fabs (x [j]) < EIG_TOL) 
	    x [j] = 0;
	}

	for (int i=1; i<np; i++)
	  for (int j=i; j<np; j++)
	    if (fabs (x [i] * x [j]) > 1e-6) {
	      coeff [nterms] = (i==j) ? 
		(x [i] * x [j]) : 
		(2 * x [i] * x [j]);
	      ind   [nterms++] = indexQ (i-1,j-1,n);
	    }

	for (int i=1; i<np; i++) 
	  if (fabs (x [i] * x [0]) > 1e-6) {
	    coeff [nterms]   = 2 * x [i] * x [0];
	    ind   [nterms++] = i-1;
	  }

	cut -> setRow (nterms, ind, coeff);
	cut -> setLb (- x [0] * x [0]);

	cs.insert (cut);

	delete cut;
	delete ind;
	delete coeff;
      }

      free (x);
    }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  /* Separator for sdp cuts.  
   *
   * Solve the problem
   *
   * min f(x) = * (1/2)*x'*A*x + b'*x
   *
   * subject to the Euclidean norm constraint
   *
   * norm (x) <= delta.
   *
   * Call the procedure dgqt in the nmtr package (thanks to Jorge
   * More' for the info)
   */

  static 
  double *A = NULL, /* Coefficient matrix of the quadratic term       */
         *b,        /* coefficient vector of the linear term          */
         *x,        /* optimal solution                               */
         *z,        /* (used by dgqt)                                 */
         *wa1,      /* (used by dgqt)                                 */
         *wa2,      /* (used by dgqt)                                 */
          mu = 0;   /* estimate (before) and real value (after) of the
		       Lagrangian multiplier of the norm constraint.
		       Static so that we set it to 0 first and keep it
		       between iterations */

  A   = (double *) malloc (np*np * sizeof (double));
  b   = (double *) malloc (np    * sizeof (double));
  x   = (double *) malloc (np    * sizeof (double));
  z   = (double *) malloc (np    * sizeof (double));
  wa1 = (double *) malloc (np    * sizeof (double));
  wa2 = (double *) malloc (np    * sizeof (double));

  register int i,j;

  /* set to zero all entries of A */
  for (j=np, i=j*=j; i--;) 
    *A++ = 0.;
  A -= j;

  /* and those of b */
  for (i=j=np; i--;) 
    *b++ = 0.;
  b -= j;

  int    niter = 0;
  double f     = 0;

  while (niter < MAX_OUT_ITER) {

    /* fill in matrix A and vector b */
    {
      register int i,j;

      A [0] = 1;

      for (j=1; j<np; j++)
	A [np*j] = sol [j-1];

      b [0] = 0;

      for (i=1; i<np; i++) {

	double bo = drand48 ();
	/*	printf ("::: %.4f <=> %.4f\n", bo, 1. - (double) niter / MAX_OUT_ITER);*/

	if (bo <= (1. - (double) niter / MAX_OUT_ITER))
	  for (j=i; j<np; j++)
	    //A [n*j+i] = Q [i] [j];
	    A [np*j+i] = sol [indexQ (i-1,j-1,n)];
	else {
	  /*	  printf ("::: %d skip %d\n", niter, i);*/
	  j=i;
	  A [np*(j++) + i] = 1;
	  for (; j<n; j++)
	    A [np*j+i] = 0;
	}

	//	b [i] = 0;//- sol [i];
      }
    }

    double 
      RTOL = 1e-5,
      ATOL = 1e-5,
      rhs  = 1;

    int lda,          /* "leading" dimension (as mentioned in dgqt.f) */
      dim = lda = np, /* dimension of vector x                        */ 
      NITER = 1000,   /* # iterations                                 */
      info = 3;       /* output status of the procedure (see dgqt.f)  */

    /****************************** call the TRM procedure ************************/
    dgqt_ (&dim,                 /* dimension of x                                */
	   A, 		         /* quadratic coefficient matrix                  */
	   &lda, 	         /* leading dimension of A                        */
	   b, 		         /* linear coefficient vector                     */
	   &rhs,   	         /* right-hand side of norm constraint            */
	   &RTOL, &ATOL, &NITER, /* rel./abs. tolerance, max iterations           */
	   &mu, 	         /* multiplier of norm constraint                 */
	   &f,                   /* optimum of the problem                        */ 
	   x, 		         /* optimal solution                              */
	   &info, 	         /* exit status                                   */
	   z, wa1, wa2);         /* vectors used by dgqt                          */
    /******************************************************************************/


#define DO_DEBUG 0

    if (DO_DEBUG) {
      int i;
      printf ("solution: %.4f; x:\n", f);
      for (i=0; i<np; i++) 
	printf ("%+.3f ", x [i]);
      printf ("\n");
    }

    switch (info) {
    case 1:  if (DO_DEBUG) printf (" :: dgqt: OK, RELATIVE accuracy\n");   break;
    case 2:  if (DO_DEBUG) printf (" :: dgqt: OK, ABSOLUTE accuracy\n");   break;
    case 3:  printf (" :: dgqt: precision issues\n");        break;
    case 4:  printf (" :: dgqt: error, cannot terminate\n"); break;
    default: printf (" :: dgqt: undefined return status\n");
    }

    niter++;

    if (info <= 2) break;
  }

  // create corresponding cut ////////////////////////////////////

  if (f < 0) { // there is a violated cut

    OsiRowCut *cut = new OsiRowCut;

    double *coeff  = new double [N];
    int    *ind    = new int    [N];
    int     nterms = 0;

    for (int i=1; i<np; i++)
      for (int j=i; j<np; j++)
	if (fabs (x [i] * x [j]) > 1e-6) {
	  coeff [nterms] = (i==j) ? 
	    (x [i] * x [j]) : 
	    (2 * x [i] * x [j]);
	  ind   [nterms++] = indexQ (i-1,j-1,n);
	}

    for (int i=1; i<np; i++) 
      if (fabs (x [i] * x [0]) > 1e-6) {
	coeff [nterms]   = 2 * x [i] * x [0];
	ind   [nterms++] = i-1;
      }

    cut -> setRow (nterms, ind, coeff);
    cut -> setLb (- x [0] * x [0]);
    cs.insert (cut);

    delete cut;
    delete ind;
    delete coeff;
  }

  // cleanup /////////////////////////////////////////////////////

  free (A); free (b);
  free (z); free (x);
  free (wa1); free (wa2);

  ncalls++;
  return;
}
