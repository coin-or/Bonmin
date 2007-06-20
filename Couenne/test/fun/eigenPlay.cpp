#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>


void SdpCutGen::eigenPlay (OsiCuts &cs, int n, int m, double *vector, double *value) const {

  double *x = (double *) malloc (n * sizeof (double));

  for (int k=0; k<m; k++) { // eigenvalue negative enough

    /*
      printf ("%3d: %6.3f -- ", k, w [k]);
      for (int j=0; j<np; j++) 
      printf ("%.3f ", z [k*np + j]);
      printf ("\n");
    */

    double *zbase = vector + k * n;

#define EIG_TOL 0

    for (int j=0; j<n; j++) {
      x [j] = *zbase++;
      if (fabs (x [j]) < EIG_TOL) 
	x [j] = 0;
    }

    genSDPcut (cs, x, x);
  }

  free (x);
}
