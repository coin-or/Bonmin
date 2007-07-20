/*
 * Name:    exprBQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: method to compute value of an expr?BQuad
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#include <exprBQuad.hpp>

CouNumber computeQBound (int sign, exprQuad *e) {

  return (sign < 0) ? -COUENNE_INFINITY : COUENNE_INFINITY;

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
    *ub = expression::Lbounds (),
    bound = e -> getc0 ();

  if (sign < 0) {

    while (nlt--) {

      CouNumber coe = *lc++;
      int       ind = *li++;

      if (coe < 0) bound += coe * lb [ind];
      else         bound += coe * ub [ind];
    }

    while (nqt--) {

      int       i = *qi++,
                j = *qj++;

      CouNumber coe = *qc++,
  	        lbi = lb [i],
  	        ubi = ub [i],
  	        lbj = lb [j],
         	ubj = ub [j],
 	        b1 = coe * lbi * lbj,
 	        b2 = coe * lbi * ubj,
 	        b3 = coe * ubi * lbj,
 	        b4 = coe * ubi * ubj;

      bound += mymin (mymin (b1, b2), mymin (b3, b4));
    }
  } else {

    while (nlt--) {

      CouNumber coe = *lc++;
      int       ind = *li++;

      if (coe > 0) bound += coe * lb [ind];
      else         bound += coe * ub [ind];
    }

    while (nqt--) {

      int       i = *qi++,
                j = *qj++;

      CouNumber coe = *qc++,
  	        lbi = lb [i],
  	        ubi = ub [i],
  	        lbj = lb [j],
         	ubj = ub [j],
 	        b1 = coe * lbi * lbj,
 	        b2 = coe * lbi * ubj,
 	        b3 = coe * ubi * lbj,
 	        b4 = coe * ubi * ubj;

      bound += mymax (mymax (b1, b2), mymax (b3, b4));
    }
  }

  return bound;

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

  /*
  int nqt = e -> getnQTerms ();

  CouNumber 
    *hL = new CouNumber [nqt],
    *hU = new CouNumber [nqt],
    *vL = new CouNumber [nqt],
    *vU = new CouNumber [nqt],
    *xl = expression::Lbounds (),
    *xu = expression::Ubounds ();

  for (register int i=nqt; i--;) 
    *hL++ = *vL++ = *hU++ = *vU++ = 0;

  hL -= nqt;  hU -= nqt;
  vL -= nqt;  vU -= nqt;
  */

  //  for (register int i = nqt; i--;) 

}
