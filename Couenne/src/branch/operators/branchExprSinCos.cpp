/*
 * Name:    branchExprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for sines/cosines
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <exprSin.hpp>
#include <exprCos.hpp>
#include <CouennePrecisions.hpp>
#include <CouenneTypes.hpp>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>


/// generalized procedure for both sine and cosine
CouNumber trigSelBranch (expression *w, 
			 const OsiBranchingInformation *info,
			 int &ind, 
			 double * &brpts, 
			 int &way,
			 enum cou_trig type) {

  // for now, apply default branching rule

  ind = -1;
  return 0.;
}
