/*
 * Name:    exprBQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: method to compute value of an expr?BQuad
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <exprBQuad.hpp>

//#define DEBUG

CouNumber computeQBound (int sign, exprQuad *e) {

  //  return (sign > 0) ? COUENNE_INFINITY : -COUENNE_INFINITY;

  // compute lower (if sign == -1) or upper (sign == +1) bound of an
  // exprQuad based on the information obtained through
  // alpha-convexification, if any, or as follows:
  //
  // w = a0 + a'x + x'Qx = 
  //   = a0 + sum{i} [(a_i + sum {j} q_ij * x_j) * x_i] =
  //   = a0 + sum{i} [                       z_i * x_i]
  // 
  // So some bound on z_i can be computed and a bound on the whole
  // expression should be better than what can be obtained by summing
  // all bounds separately.
  //
  // Notice that the above computation is fast and may be better than
  // the convexification after some updates in the variable bounds
  // without updating the convexification. Notice also that the
  // direction can also be vertical, not only horizontal

  int 
    nlt = e -> getnLTerms (),
    *li = e -> getIndices (),

    nqt = e -> getnQTerms (),
    *qi = e -> getQIndexI (),
    *qj = e -> getQIndexJ ();

  CouNumber
    *lc = e -> getCoeffs  (),
    *qc = e -> getQCoeffs (),
    *lb = expression::Lbounds (),
    *ub = expression::Ubounds (),
    bound = e -> getc0 (),
    term;

#ifdef DEBUG
  printf ("\n");
  for (int i=0; i<12; i++) printf ("%3d [%g,%g]\n",i, lb [i], ub [i]);
  e -> print ();
  printf ("\n (%g)\n ", bound);
#endif

  if (sign < 0) { // compute lower bound ////////////////////////////////////////////////

    while (nlt--) {

#ifdef DEBUG
      printf ("lin %d %g %g ", *li, *lc, (*lc < 0) ? ub [*li] : lb [*li]);
#endif

      if (*lc < 0) {if ((term = ub [*li++]) >  COUENNE_INFINITY) return -COUENNE_INFINITY;} 
      else         {if ((term = lb [*li++]) < -COUENNE_INFINITY) return -COUENNE_INFINITY;}

      bound += *lc++ * term;

#ifdef DEBUG
      printf (" --> %g\n", bound);
#endif
    }

    while (nqt--) {

      int i = *qi++,
          j = *qj++;

      CouNumber 
	coe = *qc++,
	lbi = lb [i],
	ubi = ub [i];

      if (i==j) {

	if (coe > 0) term = (ubi <= 0) ? (ubi * ubi) : (lbi >= 0) ? (lbi * lbi) : 0;
	else if ((term = mymax (lbi*lbi, ubi*ubi)) > COUENNE_INFINITY) 
	  return -COUENNE_INFINITY;

	term *= coe;

#ifdef DEBUG
	printf ("Qii %d %g %g -> %g\n", i, coe, term, bound + term);
#endif
      } else {

	coe *= 2;

	CouNumber
	  lbj = lb [j], ubj = ub [j],
	  b1 = coe * lbi * lbj, 
	  b2 = coe * lbi * ubj,
	  b3 = coe * ubi * lbj, 
	  b4 = coe * ubi * ubj;

	if (fabs (lbi) == 0) b1 = b2 = 0;
	if (fabs (lbj) == 0) b1 = b3 = 0;
	if (fabs (ubi) == 0) b3 = b4 = 0;
	if (fabs (ubj) == 0) b2 = b4 = 0;

	if ((term = mymin (mymin (b1, b2), mymin (b3, b4))) < -COUENNE_INFINITY) 
	  return -COUENNE_INFINITY; 

#ifdef DEBUG
	printf ("Qij %d %d %g %g -> %g\n", i, j, coe, term, bound + term);
#endif
      }

      //      if ((i!=j) || (lbi >= 0) || (ubi <= 0))
      bound += term;
    }
  } else { // compute upper bound /////////////////////////////////////////////////////////////

    while (nlt--) { // linear part

#ifdef DEBUG
      printf ("lin %d %g %g %g\n", *li, *lc, (*lc < 0) ? ub [*li] : lb [*li], bound);
#endif

      if (*lc > 0) {if ((term = ub [*li++]) >  COUENNE_INFINITY) return COUENNE_INFINITY;}
      else         {if ((term = lb [*li++]) < -COUENNE_INFINITY) return COUENNE_INFINITY;}

      bound += *lc++ * term;
    }

    while (nqt--) { // quadratic part

      int i = *qi++,
          j = *qj++;

      CouNumber 
	coe = *qc++,
	lbi = lb [i], 
	ubi = ub [i];

      if (i==j) {

	if (coe > 0) term = mymax (lbi * lbi, ubi * ubi);
	else         term = (ubi <= 0) ? (ubi * ubi) : (lbi >= 0) ? (lbi * lbi) : 0;

	if (term > COUENNE_INFINITY) 
	  return COUENNE_INFINITY;

	term *= coe;

#ifdef DEBUG
	printf ("Qii %d %g %g --> %g ", i, coe, term, bound + term);
#endif
      } else {

	coe *= 2;

	CouNumber 
	  lbj = lb [j], 
	  ubj = ub [j],
	  b1 = coe * lbi * lbj,
	  b2 = coe * lbi * ubj,
	  b3 = coe * ubi * lbj,
	  b4 = coe * ubi * ubj;

	// I hate this... but when you see 
	//
	//   mymax (mymax (-0, nan), mymax (nan, -inf)) =
	// = mymax (nan, -inf) = -inf
	//
	// you feel changed.

	if (fabs (lbi) == 0) b1 = b2 = 0;
	if (fabs (lbj) == 0) b1 = b3 = 0;
	if (fabs (ubi) == 0) b3 = b4 = 0;
	if (fabs (ubj) == 0) b2 = b4 = 0;

	if ((term = mymax (mymax (b1, b2), mymax (b3, b4))) > COUENNE_INFINITY)
	  return COUENNE_INFINITY;

#ifdef DEBUG
	printf ("Qij %d %d %g %g -> %g [%g,%g,%g,%g] (%g,%g,%g,%g) {%g,%g,%g}\n", 
		i, j, coe, term, bound + term,
		lbi, ubi, lbj, ubj, 
		b1, b2, b3, b4,
		mymax (b1, b2), mymax (b3, b4),
		mymax (mymax (b1, b2), mymax (b3, b4)));
#endif
      }

      bound += term;
    }
  }

  return bound;
}
