/*
 * Name:    conv-exprPow-envelope.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.hpp>
#include <rootQ.hpp>
#include <exprPow.hpp>
#include <CouennePrecisions.hpp>
#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

#define COUENNE_POW_CONST 1
#define POW_FACTOR 2

// A unary function to be given to addEnvelope with a local exponent,
// given here as a static variable
static CouNumber exponent;

// function x^k
inline static CouNumber power_k (CouNumber x) 
{return safe_pow (x, exponent);}


// function k*x^(k-1)
inline static CouNumber power_k_prime (CouNumber x) 
{return exponent * safe_pow (x, exponent-1);}


// function k*(k-1)*x^(k-2)
inline static CouNumber power_k_dblprime (CouNumber x) 
{return exponent * (exponent-1) * pow (x, exponent-2);}


// adds convex (upper/lower) envelope to a power function

void addPowEnvelope (const CouenneCutGenerator *cg, OsiCuts &cs,
		     int wi, int xi,
		     CouNumber x, CouNumber y,
		     CouNumber k, 
		     CouNumber l, CouNumber u,
		     int sign) {
  exponent = k;

  // set x to get a deeper cut (so that we get a tangent which is
  // orthogonal with line through current- and tangent point)

  if (!(cg -> isFirst ()))
    x = powNewton (x, y, power_k, power_k_prime, power_k_dblprime);

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

  // convex envelope
  cg -> addEnvelope (cs, sign, power_k, power_k_prime, 
		     wi, xi, x, l, u);
}
