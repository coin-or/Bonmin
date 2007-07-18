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

  //  for (register int i = nqt; i--;) 

}
