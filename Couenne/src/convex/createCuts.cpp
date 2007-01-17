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

OsiRowCut *CouenneCutGenerator::createCut (CouNumber rhs, int sign, 
					   int i1, CouNumber c1,
					   int i2, CouNumber c2,
					   int i3, CouNumber c3) const {

  // a maximum of three terms are allowed here. Index -1 means the
  // term is not considered

  int nterms = 0;

  if (i1 >= 0) nterms++; // useless, but you never know...
  if (i2 >= 0) nterms++;
  if (i3 >= 0) nterms++;

  bool check = isFirst () || !(addViolated ());

  if (!check) { // need to check violation 

    // compute violation

    CouNumber violation = - rhs + c1 * X (i1);
    if (i2 >= 0) violation     += c2 * X (i2);
    if (i3 >= 0) violation     += c3 * X (i3);

    if      ((violation >   COUENNE_EPS) && (sign <= 0)) check = true;
    else if ((violation < - COUENNE_EPS) && (sign >= 0)) check = true;
  }

  if (check) { // that is, if this is the first call, if we also want
	       // unviolated cuts, or if it is violated

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

    return cut;
  }
  else return NULL;
}
