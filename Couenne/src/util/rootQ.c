/*
 * Name:    rootQ.c
 * Author:  Pietro Belotti
 * Purpose: find roots of polynomial Q^k(x) (see Liberti and Pantelides, 2003)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.h>
#include <CouennePrecisions.h>
#include <stdio.h>

/*
inline int round (double x) {
  return (int) (ceil (x - 0.5));
}
*/

/* compute Q(x)*/

CouNumber Q (int k, CouNumber x) {

  register CouNumber xp;
  register int i;
  register CouNumber Q;

  k *= 2;

  xp = x;
  Q = 1;

  for (i=2; i<=k; i++) {

    Q += (CouNumber) i * xp;
    xp *= x;
  }

  return Q;
}


/*
 * Find roots of polynomial $Q^k(x) = \sum_{i=1}^{2k} i x^{i-1}$. Used
 *  in convexification of powers with odd exponent
 */

CouNumber rootQ (int k) {

  if (k==1) return - 0.5;
  else {

    register CouNumber l = - 1.0 + 0.5 / k, 
                    u = - 0.5,
                    Ql = Q (k, l), Qu = Q (k, u), Qm,
                    midpoint;
    do {

      midpoint = 0.5 * (l+u);/* (- Ql * u + Qu * l) / (Qu - Ql); */
      Qm = Q (k, midpoint);

      /*      printf ("[%.4f, %.4f] --> %.4f: %.24f\n", l, u, midpoint, Qm); */

      if (Qm<0) {l = midpoint; Ql = Qm;}
      else      {u = midpoint; Qu = Qm;}

    } while (fabs (Qm) > 1e-15);

    return midpoint;
  }
}

#ifdef DEBUG_ROOTQ
int main () {

  register int k;
  CouNumber x, q;

  for (k=30; --k;) {

    printf ("root, %3d -> %.15f\n", 2*k+1, rootQ (k));
    /*
    printf ("k=%3d: ", 2*k+1);
    for (x = -1.0; x < 0.4; x += 0.1) {
      //      Q (k, x, &q, NULL, NULL); 
      printf ("[%.2f, %.3f] ", x, rootQ);
    }
    printf ("\n");
    */
  }
}
#endif
