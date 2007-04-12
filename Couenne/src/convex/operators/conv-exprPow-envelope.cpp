/*
 * Name:    conv-exprPow-envelope.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.h>
#include <rootQ.h>
#include <exprPow.h>
#include <CouennePrecisions.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

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


// adds convex (upper/lower) envelope to a power function

void addPowEnvelope (const CouenneCutGenerator *cg, OsiCuts &cs,
		     int wi, int xi,
		     CouNumber x, CouNumber k, 
		     CouNumber l, CouNumber u,
		     int sign) {

  int ns = cg -> nSamples ();

  exponent = k;
  /*
  if (k <= 0) k = 1;

  if (l < - COUENNE_INFINITY + 1) {
    if (u > COUENNE_INFINITY - 1) {

      l = - pow (POW_FACTOR, 1. / (1. + k)) * ns/2;
      u =   pow (POW_FACTOR, 1. / (1. + k)) * ns/2;
    } else 
      l = u - pow (POW_FACTOR, 1. / (1. + k)) * ns;

  } else if (u > COUENNE_INFINITY - 1) 
    u = l + pow (POW_FACTOR,   1. / (1. + k)) * ns;
  */

  // limit the bounds for the envelope

  if (l < - COUENNE_INFINITY + 1) {
    if (u > COUENNE_INFINITY - 1) {
      l = x - 1;
      u = x + 1;
    } else 
      l = x - 1;
  } else 
    if (u > COUENNE_INFINITY - 1) 
      u = x + 1;

  // convex envelope

  cg -> addEnvelope (cs, sign, power_k, power_k_prime, 
		     wi, xi, x, l, u, cg -> isFirst ());
}
