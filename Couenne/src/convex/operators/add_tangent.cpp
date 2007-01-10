/*
 * Name:    add_tangent.C
 * Author:  Pietro Belotti
 * Purpose: a method to add an OsiRowCut of the form w >= ax + b
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprAux.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// add half-plane defined by x, w, and slope, and <= or >= according to sign

void CouenneCutGenerator::addTangent (OsiCuts &cs, int wi, int xi, 
		 CouNumber x, CouNumber w, CouNumber slope, int sign) const { 


  CouNumber violation = X (wi) - w - slope * (X (xi) - x);
  bool violated = 
    ((violation >   COUENNE_EPS) && sign <= 0) ||
    ((violation < - COUENNE_EPS) && sign >= 0);

  // add tangent only if first call or if violated

  if (firstcall_ || !addviolated_ || violated) {

    CouNumber *coeff = new CouNumber [2];
    int       *index = new int       [2];
    OsiRowCut *cut   = new OsiRowCut;

    coeff [0] = 1.;      index [0] = wi;
    coeff [1] = - slope; index [1] = xi;

    if (sign >= 0) cut -> setLb (w - slope * x);
    if (sign <= 0) cut -> setUb (w - slope * x);

    cut -> setRow (2, index, coeff);

    printf ("Tangent: "); cut -> print ();

    cs.insert (cut);
  }
}


// add half-plane defined by two points (x1,y1) and (x2, y2), and
// sign. Sign is only valid if x1 < x2

