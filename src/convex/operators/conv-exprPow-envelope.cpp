/*
 * Name:    conv-exprPow-envelope.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "CouenneTypes.hpp"
#include "rootQ.hpp"
#include "exprPow.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"
#include "funtriplets.hpp"


// adds convex (upper/lower) envelope to a power function

void addPowEnvelope (const CouenneCutGenerator *cg, OsiCuts &cs,
		     int wi, int xi,
		     CouNumber x, CouNumber y,
		     CouNumber k, 
		     CouNumber l, CouNumber u,
		     int sign) {

  // set x to get a deeper cut (so that we get a tangent which is
  // orthogonal with line through current- and tangent point)

  if (!(cg -> isFirst ())) {

    powertriplet pt (k);
    x = powNewton (x, y, &pt);
  }

  if      (x<l) x=l;
  else if (x>u) x=u;

  // limit the bounds for the envelope

  CouNumber powThres = (k<=1) ? COU_MAX_COEFF: pow (COU_MAX_COEFF, 1./k),
            step     = (1 + log (1. + (double) (cg -> nSamples ()))) * powThres / COU_MAX_COEFF;

  if (l < - powThres + 1) {
    l = x - step;
    if (u > powThres - 1)
      u = x + step;
  } else 
    if (u > powThres - 1) 
      u = x + step;

  powertriplet pt (k);

  // convex envelope
  cg -> addEnvelope (cs, sign, &pt, 
		     wi, xi, x, l, u);
}
