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
#include <CouenneCutGenerator.h>


// general procedure for inserting a linear cut with up to three
// variables

OsiRowCut *CouenneCutGenerator::createCut (CouNumber rhs, int sign, 
					   int i1, CouNumber c1,
					   int i2, CouNumber c2,
					   int i3, CouNumber c3,
					   bool is_global) const {

  // a maximum of three terms are allowed here. Index -1 means the
  // term is not considered

  int nterms = 0;

  if (i1 >= 0) nterms++; // useless, but you never know...
  if (i2 >= 0) nterms++;
  if (i3 >= 0) nterms++;

  if (!firstcall_ && addviolated_) { // need to check violation 

    // compute violation

    CouNumber violation = - rhs + c1 * X (i1);
    if (i2 >= 0) violation     += c2 * X (i2);
    if (i3 >= 0) violation     += c3 * X (i3);

    // return NULL if not violated
    if (((violation <   COUENNE_EPS) || (sign > 0)) &&
	((violation > - COUENNE_EPS) || (sign < 0)))
      return NULL;
  }

  // You are here if:
  //
  // 1) this is the first call to CouenneCutGenerator::generateCuts()
  // 2) we also want unviolated cuts
  // 3) the cut is violated

  CouNumber *coeff = new CouNumber [nterms]; 
  int       *index = new int       [nterms];
  OsiRowCut *cut   = new OsiRowCut;

  coeff [0] = c1; index [0] = i1;
  if (i2 >= 0) {coeff [1] = c2; index [1] = i2;}
  if (i3 >= 0) {coeff [2] = c3; index [2] = i3;}

  if (sign <= 0) cut -> setUb (rhs);
  if (sign >= 0) cut -> setLb (rhs);

  cut -> setRow (nterms, index, coeff);

  delete [] coeff;
  delete [] index;

  // some convexification cuts (as the lower envelopes of convex
  // functions) are global, hence here is a tool to make them valid
  // throughout the code

  cut -> setGloballyValid (is_global);

  return cut;
}
