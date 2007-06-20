/*
 * Name:    partialLU.c
 * Author:  Pietro Belotti
 * Purpose: decompose matrix with nonzero pivots on the diagonal and
 * leave all null pivots in principal minor on the lower right
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License.
 */

#include <stdio.h>
#include <stdlib.h>

#include <CouennePrecisions.h>

/*
 *  Check if matrix is positive semidefinite by eliminating
 *  non-diagonal elements. Use only pivots on the diagonal and permute
 *  both rows and columns
 *
 *
 *
 *
 *
 */

double **partialLU (double **Q, int n) {

  /* While at least one element on the diagonal is nonzero */
  /*   Use that element as pivot to eliminate column       */
  /*   store coefficient in L                              */
  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                                                       */
  /*                                                       */

  int i, *perm = (int *) malloc (n * sizeof (int));

  /* initialize permutation */
  for (i=0; i<n;) *perm++ = i++;
  perm -= n;

  for (i=0; i<n; i++) {

    register int j=i;

    for (; j<n;) { /* look for nonzero element on the remaining diagonal */
      int p = perm [j++];
      if (fabs (Q [p] [p]) > COUENNE_EPS) 
	break;
    }

    if (j==n) break;

    if (j>i) { /* had to go down a little to find nonzero, hence permute */
      int swap = perm [j];
      perm [j] = perm [i];
      perm [i] = swap;
    }

    /* eliminate everything under element i */


  }

  return NULL;
}
