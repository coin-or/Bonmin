#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

/*****************************************/
int main (int argc, char **argv) {


  struct timeval tv;
  gettimeofday (&tv, NULL);

  /*  srand48 (tv.tv_sec * 1e6 + tv.tv_usec);*/
  srand48 (tv.tv_usec);

  int n = atoi (argv [1]), i, j;
  double *b, **Q, 
    l = atof (argv [2]),
    u = atof (argv [3]);

  printf ("%d\n", n);

  b = (double *)  malloc (n * sizeof (double));
  Q = (double **) malloc (n * sizeof (double *));

  for (i=0; i<n; i++)
    Q [i] = (double *) malloc (n * sizeof (double));

  for (j=0; j<n; j++)
    printf ("%5.0f ", l + drand48 () * (u-l));

  printf ("\n");

  for   (i=0; i<n; i++)
    for (j=i; j<n; j++)
      Q [i] [j] = Q [j] [i] = l + drand48 () * (u-l);

  for   (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      printf ("%5.0f ", Q [i] [j]);
    printf ("\n");
  }

  return 0;
}
