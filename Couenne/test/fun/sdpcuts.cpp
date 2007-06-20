/*
 * Apply semidefinite cuts to solve Box constrained Quadratic
 * Programming problems
 *
 * (C) Pietro Belotti, Carnegie Mellon University
 */

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <OsiCpxSolverInterface.hpp>
#include <SdpCutGen.hpp>

void populateProblem (int, double *, double **, OsiSolverInterface *);

/*****************************************/
int main (int argc, char **argv) {


  struct timeval tv;
  gettimeofday (&tv, NULL);

  /*  srand48 (tv.tv_sec * 1e6 + tv.tv_usec);*/
  srand48 (tv.tv_usec);

  int n, i, j;
  double *b, **Q;
  FILE *f = fopen (argv [1], "r");

  /* read instance ****************************************************/

  if (!fscanf (f, "%d", &n) || (n <= 0)) return -1;

  b = (double *)  malloc (n * sizeof (double));
  Q = (double **) malloc (n * sizeof (double *));

  for (i=0; i<n; i++)
    Q [i] = (double *) malloc (n * sizeof (double));

  for (j=0; j<n; j++)
    fscanf (f, "%lf", b + j);

  for   (i=0; i<n; i++)
    for (j=0; j<n; j++)
      fscanf (f, "%lf", Q [i] + j);

  fclose (f);

  /***** start solver ***********************************************/

  OsiCpxSolverInterface si;

  populateProblem (n, b, Q, &si);

  SdpCutGen scg (n, b, Q);

  si.initialSolve();

  int niter = 0, ncuts, ntotcuts=0;

  do {

    struct rusage use;
    double time;

    getrusage (RUSAGE_SELF, &use);
    time = use.ru_utime.tv_sec + 1e-6 * use.ru_utime.tv_usec;

    printf ("::: %4d %8.2f %5d %12.4f\n", 
	    niter, time, ntotcuts, si. getObjValue());

    OsiCuts cs;
    scg.generateCuts (si, cs);

    si.applyCuts (cs);
    si.resolve();

    niter++;

    ntotcuts += ncuts = cs.sizeRowCuts ();

    //    char fname [20];
    //    sprintf (fname, "sdp_%d", niter);
    //    si.writeLp (fname);

  } while ((ncuts > 0) && (niter < 500));

  //  si.writeLp ("sdp-cut");

  /* cleanup */

  free (b);
  for (i=0; i<n; i++)
    free (Q [i]);
  free (Q);

  return 0;
}
