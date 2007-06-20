#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

// sdpcut separator
void SdpCutGen::separateBK (const OsiSolverInterface &si, 
			    OsiCuts &cs) const {

  int n  = *(const_cast <int *> (&n_)), 
      np = n + 1;

  const double *sol = si.getColSolution ();

  /* 
   *  get decomposition of 
   *
   *          / 1  x^T \
   *  X'  =   \ x  X   /
   *
   *  with both row and column operations using ONLY pivots on the
   *  diagonal
   */

  double **X = (double **) malloc (np * sizeof (double *));

  for (int i=0; i<np; i++) {

    register int j=i;

    X [i] = (double *) malloc (np * sizeof (double));

    X [i] [i] = sol [indexQ (i-1,i-1,n)];

    for (++j; j<np; j++) 
      X [i] [j] = 
      X [j] [i] = sol [indexQ (i-1,j-1,n)];
  }

  /*
  register int i,j;

    A [0] = 1;

    for (j=1; j<np; j++)
      A [np*j] = sol [j-1];

    for (i=1; i<np; i++)
      for (j=i; j<np; j++)
	A [np*j+i] = sol [indexQ (i-1,j-1,n)];
  }
  */

  double **list = partialLU (X, np);

  for (int i=0; list [i]; i++)
    genSDPcut (cs, list [i], list [i]);

  /*
  nrows = 3;
  ncols = 3;
  A [0] =   0; A [3] = 9; A [6] =  0;
  A [1] =   0; A [4] = 0; A [7] =  4;
  A [2] =  -5; A [5] = 0; A [8] =  0;
  lda = 3;
  */
}
