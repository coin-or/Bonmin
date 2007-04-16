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


// general procedure for inserting a linear cut with up to three
// variables. Return 1 if cut inserted, 0 if none, <0 if error

int CouenneCutGenerator::createCut (OsiCuts &cs,
				    CouNumber rhs, int sign, 
				    int i1, CouNumber c1,
				    int i2, CouNumber c2,
				    int i3, CouNumber c3,
				    bool is_global)       const {

  // a maximum of three terms are allowed here. Index -1 means the
  // term is not considered

  int nterms = 0;

  if (i1 >= 0) nterms++; // useless, but you never know...
  if (i2 >= 0) nterms++;
  if (i3 >= 0) nterms++;

  if (!nterms) // nonsense cut
    return 0;

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

    OsiColCut *cut = new OsiColCut;

    if (sign <= 0) cut -> setUbs (1, &i1, &rhs);
    if (sign >= 0) cut -> setLbs (1, &i1, &rhs);

    cut -> setGloballyValid (is_global); // global?

    cs.insert (cut);

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
  }

  return 1;
}
