/*
 * Name:    exprBQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: method to compute value of an expr?BQuad
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "exprBQuad.hpp"

//#define DEBUG

CouNumber exprQuad::computeQBound (int sign) {

  // 1) loose: disaggregated bound
  // 2) tighter: aggregate coefficient per variable
  // 3) tightest: solve (convex) QP on available alpha-convexification

  // 1) loose bound
  //
  // w = a0 + a'x + x'Qx means that its lower bound is
  //
  // w_l = a0 + a'_+ x_l + a'_- x_u + \sum_


  // 2) tighter bound -- TODO
  //
  // Compute lower (if sign == -1) or upper (sign == +1) bound of an
  // exprQuad based on the information obtained through
  // alpha-convexification, if any, or as follows:
  //
  // w = a0 + a'x + x'Qx = 
  //   = a0 + sum{i} [(a_i + sum {j>=i} q_ij * x_j) * x_i] =
  //   = a0 + sum{i} [                          z_i * x_i]
  // 
  // Thus, some bound on z_i can be computed and a bound on the whole
  // expression should be better than what can be obtained by summing
  // all bounds separately.
  //
  // Notice that the above computation is fast and may be better than
  // the convexification after some updates in the variable bounds
  // without updating the convexification. Notice also that the
  // direction can also be vertical, not only horizontal.

  CouNumber bound = c0_, term;

  // derive linear part (obtain constant)
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    CouNumber 
      coe = el -> second, term = 0.,
      lb  = el -> first -> lb (),
      ub  = el -> first -> ub ();

    if ((coe < 0.) && (sign < 0) || 
	(coe > 0.) && (sign > 0)) 
      {if    ((term=ub) >  COUENNE_INFINITY) return (sign < 0)? -COUENNE_INFINITY : COUENNE_INFINITY;}
    else {if ((term=lb) < -COUENNE_INFINITY) return (sign < 0)? -COUENNE_INFINITY : COUENNE_INFINITY;}

    bound += coe * term;
  }

#ifdef DEBUG
  printf ("quadBound --- linear, %cb = %g\n", (sign < 0) ? 'l' : 'u', bound);
#endif

  // derive quadratic part (obtain linear part)
  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    int xind = row -> first -> Index ();

      CouNumber 
	lbi = row -> first -> lb (),
	ubi = row -> first -> ub ();

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      int yind = col -> first -> Index ();

      CouNumber coe = col -> second;

      if (xind == yind) { // term of the form q_ii x_i^2

	if ((coe > 0.) && (sign < 0) ||
	    (coe < 0.) && (sign > 0)) 
	  term = (ubi < 0) ? (ubi * ubi) : (lbi > 0) ? (lbi * lbi) : 0.; //min{xi^2: xi in [lbi,ubi]
	else 
	  if ((term = CoinMax (lbi*lbi, ubi*ubi)) > COUENNE_INFINITY) 
	    return (sign < 0) ? -COUENNE_INFINITY : COUENNE_INFINITY;

	term *= coe;

#ifdef DEBUG
	printf ("Qii %d %g %g -> %g\n", xind, coe, term, bound + term);
#endif
      } else {

	coe *= 2;

	CouNumber
	  lbj = col -> first -> lb (),
	  ubj = col -> first -> ub (),
	  b1 = coe * lbi * lbj, 
	  b2 = coe * lbi * ubj,
	  b3 = coe * ubi * lbj, 
	  b4 = coe * ubi * ubj;

	if (fabs (lbi) == 0) b1 = b2 = 0;
	if (fabs (lbj) == 0) b1 = b3 = 0;
	if (fabs (ubi) == 0) b3 = b4 = 0;
	if (fabs (ubj) == 0) b2 = b4 = 0;

	if (sign < 0) {
	  if ((term = CoinMin (CoinMin (b1, b2), CoinMin (b3, b4))) < -COUENNE_INFINITY) 
	    return -COUENNE_INFINITY; 
	} else
	  if ((term = CoinMax (CoinMax (b1, b2), CoinMax (b3, b4))) >  COUENNE_INFINITY)
	    return  COUENNE_INFINITY;

#ifdef DEBUG
	printf ("Qij %d %d %g %g -> %g\n", xind, yind, coe, term, bound + term);
#endif
      }

      //      if ((i!=j) || (lbi >= 0) || (ubi <= 0))
      bound += term;
    }
  }

  return bound;
}
