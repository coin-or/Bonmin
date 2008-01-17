/*
 * Name:    branchExprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for sines/cosines
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "exprSin.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"


/// generalized procedure for both sine and cosine
CouNumber trigSelBranch (const CouenneObject *obj, 
			 const OsiBranchingInformation *info,
			 expression *&var,
			 double * &brpts, 
			 int &way,
			 enum cou_trig type) {

  // for now, apply default branching rule

  var = NULL;
  return 0.;

  // TODO: minarea
}
