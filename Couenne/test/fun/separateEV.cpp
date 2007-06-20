/*
 * Name:    separateEV.cpp
 * Author:  Pietro Belotti
 * Purpose: find eigenvalues/vectors of X' and use them as cuts
 *
 * (C) Pietro Belotti. This file is publised under the Common Public License.
 */

#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>


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

    *A = 1;

    for (j=1; j<np; j++)
      A [np*j] = sol [j-1];

    for (i=1; i<np; i++)
      for (j=i; j<np; j++)
	A [np*j+i] = sol [indexQ (i-1,j-1,n)];
  }

  int m;
  double *w, *z;

  double alpha = 1;

  for (;;) {

    dsyevx_wrapper (    np, A,   m, w, z);
    eigenPlay      (cs, np, sol, m, w, z);

    if ((alpha *= 0.95) < 1) 
      break;

    // change A so that lower right block is 
    //
    // alpha*(x x^T) + (1-alpha) X 
    //
    // rather than just X

    register int i,j;

    *A = 1;

    for (j=1; j<np; j++)
      A [np*j] = sol [j-1];

    for (i=1; i<np; i++)
      for (j=i; j<np; j++)
	A [np*j+i] =        alpha  * sol [indexQ (i-1,j-1,n)] 
	             + (1 - alpha) * sol [i-1] * sol [j-1];
  } 

  free (A); 
  free (w);
  free (z);

  ncalls++;
}
