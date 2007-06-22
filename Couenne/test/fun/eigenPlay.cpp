#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>


void SdpCutGen::eigenPlay (OsiCuts &cs, 
			   int n, const double *x, 
			   int m, double *value, double *vector) const {

  double *v = (double *) malloc (n * sizeof (double));

  for (int k=0; k<m; k++) { // eigenvalue negative enough

    if (value [k] < 0)
      violated_ = true;

    /*
      printf ("%3d: %6.3f -- ", k, w [k]);
      for (int j=0; j<np; j++) 
      printf ("%.3f ", z [k*np + j]);
      printf ("\n");
    */

    double *zbase = vector + k * n;

#define EIG_TOL 0

    for (int j=0; j<n; j++) {
      v [j] = *zbase++;
      if (fabs (v [j]) < EIG_TOL) 
	v [j] = 0;
    }

    genSDPcut (cs, v, v);

    // add cuts with sums of eigenvectors (consider pairs whose
    // eigenvalue is below a certain threshold)

#define EVAL_THRES -1.5e5

    for (int i=0; i<n; i++)
      if (value [i] < EVAL_THRES) {

	double *v1base = vector + i*n;

	for (int j=0; j<n; j++)
	  if (value [j] < EVAL_THRES) {

	    double *v2base = vector + j*n;

	    for (int k=0; k<n; k++)
	      v [k] = *v1base++ + *v2base++;

	    genSDPcut (cs, v, v);
	  }
      }
  }

  free (v);
}
