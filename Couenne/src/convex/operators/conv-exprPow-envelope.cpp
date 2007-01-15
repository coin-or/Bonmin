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
{return pow (x, exponent);}


// function k*x^(k-1)

inline static CouNumber power_k_prime (CouNumber x) 
{return exponent * pow (x, exponent-1);}


// adds convex (upper/lower) envelope to a power function

void addPowEnvelope (const CouenneCutGenerator *cg, OsiCuts &cs,
		     int wi, int xi,
		     CouNumber x, CouNumber k, 
		     CouNumber l, CouNumber u,
		     int sign) {

  int ns = cg -> nSamples ();

  //  printf ("Power: \n"); 

  if (l < - COUENNE_INFINITY + 1) {
    if (u > COUENNE_INFINITY - 1) {

      l = - pow (POW_FACTOR, k) * ns/2;
      u =   pow (POW_FACTOR, k) * ns/2;
    } else 
      l = u - pow (POW_FACTOR, k) * ns;

  } else if (u > COUENNE_INFINITY - 1) 
    u = l + pow (POW_FACTOR, k) * ns;

  // convex envelope

  exponent = k;

  cg -> addEnvelope (cs, sign, power_k, power_k_prime, wi, xi, x, l, u);
				    
  /*
  if ((cg -> ConvType () == UNIFORM_GRID) || cg -> isFirst ()) {

    // choose sampling points. 

    // If unbounded, re-bound using a rule of thumb where each point
    // is taken every z from the finite bound, z dependent on operator
 
    // now add tangent at each sampling point

    CouNumber sample = l, step = (u-l) / ns;

    for (int i = 0; i <= ns; i++) {

      addTangent (cs, wi, xi, sample, pow (sample, k), k * pow (sample, k-1), sign);

      sample += step;
    }
  } else 
    if (cg -> ConvType () == CURRENT_ONLY)
      addTangent (cs, wi, xi, x, pow (x, k), k * pow (x, k-1), sign);

  else {

    CouNumber sample = x;

    addTangent (cs, wi, xi, x, exp (x), exp (x), +1);

    for (int i = 0; i <= ns/2; i++) {

      sample -= (x-l) / ns;
      addTangent (cs, wi, xi, sample, pow (sample, k), k * pow (sample, k-1), sign);
    }

    sample = x;

    for (int i = 0; i <= ns/2; i++) {

      sample += (u-x) / ns;
      addTangent (cs, wi, xi, sample, pow (sample, k), k * pow (sample, k-1), +1);
    }
  }
  */
}
