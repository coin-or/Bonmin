/*
 * Name:    add_tangent.C
 * Author:  Pietro Belotti
 * Purpose: a method to add an OsiRowCut of the form w >= ax + b
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprAux.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// add half-plane defined by x, w, and slope, and <= or >= according to sign

void CouenneCutGenerator::addTangent (OsiCuts &cs, int wi, int xi, 
				      CouNumber x, CouNumber w, 
				      CouNumber slope, int sign) const {

  CouNumber violation = (X (wi) - w - slope * (X (xi) - x)) / sqrt (slope*slope + 1);
  bool violated = 
    ((violation >   COUENNE_EPS) && sign <= 0) ||
    ((violation < - COUENNE_EPS) && sign >= 0);

  //  printf ("######## P=(%.8f,%.8f) xi=%d  w_(%d)=%.8f slope=%.8f sign=%d violation=%.8f\n",
  //	  x, w, xi, wi, X (wi), slope, sign, violation);

  // add tangent only if first call or if violated

  if (firstcall_ || !addviolated_ || violated) {

    CouNumber *coeff = new CouNumber [2];
    int       *index = new int       [2];
    OsiRowCut *cut   = new OsiRowCut;

    CouNumber rhs = w - slope * x;

    coeff [0] = 1.;      index [0] = wi;
    coeff [1] = - slope; index [1] = xi;

    if (sign >= 0) cut -> setLb (rhs);
    if (sign <= 0) cut -> setUb (rhs);

    cut -> setRow (2, index, coeff);

    delete [] index;
    delete [] coeff;

    //    printf ("Tangent: "); cut -> print ();

    cs.insert (cut);
  }
}
