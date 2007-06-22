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
}


#define MAX_OUT_ITER 20
#define MAX_EIGENVALUE -1e-6

// sdpcut separator
void SdpCutGen::separateTRM (const OsiSolverInterface &si, 
			      OsiCuts &cs) const {

  static int ncalls = 0;

  int n  = *(const_cast <int *> (&n_)), 
      np = n + 1;

  const double *sol = si.getColSolution ();

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

  /*
  printf ("############################################\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      printf ("%+7.2f ", sol [indexQ (i,j,n)]);
    printf ("\n");
  }
  printf ("\n");
  for (i=0; i<n; i++)
    printf ("%+7.2f ", sol [i]);
  printf ("\n############################################\n");
  */

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


#define DO_DEBUG 1

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

  printf ("---------- violation : %.4f\n", f);

  if (f < -1e-5) { // there is a violated cut
    violated_ = true;
    genSDPcut (cs, x, x);
  }

  // cleanup /////////////////////////////////////////////////////

  free (A); free (b);
  free (z); free (x);
  free (wa1); free (wa2);

  ncalls++;
}
