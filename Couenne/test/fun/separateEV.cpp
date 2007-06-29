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

  const double 
    *sol  = si.getColSolution (),
    *best = bestSol_;

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

  double *A = (double *) malloc (np*np * sizeof (double)), 
    *w, *z, 
    alpha = 0.9;

  for (int iter = 0; alpha > 0.2; alpha *= 0.3, iter++) {

    // change A so that lower right block of X is 
    //
    // alpha*(x x^T) + (1-alpha) X 
    //
    // rather than just X

    *A = 1;

    for (register int j=1; j<np; j++)
      A [np*j] = sol [j-1];

    if (best && iter)
      for (register int i=1; i<np; i++)
	for (register int j=i; j<np; j++)
	  A [np*j+i] =        alpha  * sol [indexQ (i-1, j-1, n)] 
	               + (1 - alpha) * best [i-1] * best [j-1];
    else 
      if (!iter)
	for (register int i=1; i<np; i++)
	  for (register int j=i; j<np; j++)
	    A [np*j+i] = sol [indexQ (i-1, j-1, n)];
      else break;

    int m;
    dsyevx_wrapper (    np, A,   m, w, z);
    eigenPlay      (cs, np, sol, m, w, z);

  } 

  free (A); 
  free (w);
  free (z);

  ncalls++;
}
