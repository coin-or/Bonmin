/*
 * Name:    branchExprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch suggestion for exprAbs
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "exprAbs.hpp"
#include "CouenneObject.hpp"

static const double sqrt_2 = sqrt (2.);

/// set up branching object by evaluating branching points for each
/// expression's arguments. For an exprAbs, simply branch at zero.
CouNumber exprAbs::selectBranch (const CouenneObject *obj,
				 const OsiBranchingInformation *info,
				 expression * &var,
				 double * &brpts,
				 double * &brDist, // distance of current LP
						   // point to new convexifications
				 int &way) {
  var = argument_;

  int ind = var -> Index ();

  assert ((ind >= 0) && (obj -> Reference () -> Index () >= 0));

  CouNumber x0 = info -> solution_ [ind],
            y0 = info -> solution_ [obj -> Reference () -> Index ()];

  brpts = (double *) realloc (brpts, sizeof (double));

  // the best branching point for |x| is 0, as the two subproblems
  // will have exact convexifications (lines)
  *brpts = 0.;

  way = TWO_RAND; // don't care which subtree to visit first

  // no need to compute two distances for pseudocost, as this object
  // will only branch once...

  brDist = (double *) realloc (brDist, 2 * sizeof (double));

  assert ((y0 >= x0) && (y0 >= -x0));

  brDist [0] = (x0 + y0) / sqrt_2;
  brDist [1] = (y0 - x0) / sqrt_2;

  // exact distance between current point and the two subsequent
  // convexifications
  return CoinMin (brDist [0], brDist [1]);
}
