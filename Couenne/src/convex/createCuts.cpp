/*
 * Name:   createCuts.cpp
 * Author: Pietro Belotti
 * Purpose: a standard cut creator for use with convexification
 *
 *
 * (C) 2006 Pietro Belotti, all rights reserved.
 * This file is distributed under the Common Public License
 */

#include <OsiRowCut.hpp>

#include <CouenneTypes.h>
#include <CouennePrecisions.h>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>


#define MAX_COEFF 1e10

/// general procedure for inserting a linear cut with up to three
/// variables. Return 1 if cut inserted, 0 if none, <0 if error

int CouenneCutGenerator::createCut (OsiCuts &cs,
				    CouNumber rhs, int sign, 
				    int i1, CouNumber c1,
				    int i2, CouNumber c2,
				    int i3, CouNumber c3,
				    bool is_global)       const {
  bool numerics = false;

  // a maximum of three terms are allowed here. Index -1 means the
  // term is not considered

  int nterms = 0;

  if (i1 >= 0) {if (fabs (c1) > MAX_COEFF) numerics = true; nterms++;} else c1 = 0;
  if (i2 >= 0) {if (fabs (c2) > MAX_COEFF) numerics = true; nterms++;} else c2 = 0;
  if (i3 >= 0) {if (fabs (c3) > MAX_COEFF) numerics = true; nterms++;} else c3 = 0;

  if (!nterms) // nonsense cut
    return 0;

  // cut has large coefficients/rhs, bail out
  if (numerics || (fabs (rhs) > MAX_COEFF)) {
    printf ("Warning, too large coefficients/rhs: %g, %g, %g; %g\n", c1, c2, c3, rhs);
    return 0;
  }

  if (!firstcall_ && addviolated_) { // need to check violation 

    CouNumber *x = const_cast <CouNumber *> (problem_ -> X ());

    // compute violation
    CouNumber violation = - rhs;

    if (i1 >= 0) violation += c1 * x [i1];
    if (i2 >= 0) violation += c2 * x [i2];
    if (i3 >= 0) violation += c3 * x [i3];

    // return NULL if not violated

    if (((violation <   COUENNE_EPS) || (sign > 0)) &&
	((violation > - COUENNE_EPS) || (sign < 0)))
      return 0;
  }

  // You are here if:
  //
  // 1) this is the first call to CouenneCutGenerator::generateCuts()
  // 2) we also want unviolated cuts
  // 3) the cut is violated

  // two cases: cut is of the form w1 [<|>]= alpha, hence a column
  // cut, or it is of the form (a w1 + b w2 + c w3 [<|>]= alpha), a
  // row cut

  if ((i2 < 0) && (i3 < 0)) { // column cut

    if ((fabs (c1) < COUENNE_EPS) && (fabs (rhs) > 1e10*COUENNE_EPS)) {
      printf ("nonsense column cut: %e w_%d <>= %e\n", c1, i1, rhs);
      return 0;
    }

    OsiColCut *cut = new OsiColCut;

    CouNumber bound = rhs/c1;

    if (c1 < 0) sign = -sign;

    if (sign <= 0) cut -> setUbs (1, &i1, &bound);
    if (sign >= 0) cut -> setLbs (1, &i1, &bound);

    cut -> setGloballyValid (is_global); // global?
    cs.insert (cut);
    delete cut;

  } else { // row cut

    CouNumber *coeff = new CouNumber [nterms]; 
    int       *index = new int       [nterms];
    OsiRowCut *cut   = new OsiRowCut;

    if (i1 >= 0) {coeff [0] = c1; index [0] = i1;}
    if (i2 >= 0) {coeff [1] = c2; index [1] = i2;}
    if (i3 >= 0) {coeff [2] = c3; index [2] = i3;}

    if (sign <= 0) cut -> setUb (rhs);
    if (sign >= 0) cut -> setLb (rhs);

    cut -> setRow (nterms, index, coeff);

    delete [] coeff;
    delete [] index;

    // some convexification cuts (as the lower envelopes of convex
    // functions) are global, hence here is a tool to make them valid
    // throughout the BB tree

    cut -> setGloballyValid (is_global);
    cs.insert (cut);
    delete cut;
  }

  return 1;
}
