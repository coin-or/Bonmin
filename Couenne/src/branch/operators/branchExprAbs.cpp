/*
 * Name:    branchExprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch suggestion for exprAbs
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprAbs.hpp>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneBranchingObject.hpp>
#include <CouenneObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments. For an exprAbs, simply branch at
/// zero.

CouNumber exprAbs::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {
  ind = argument_ -> Index ();

  if (ind < 0) {
    printf ("selectBranch: argument of exprAbs has negative index\n");
    exit (-1);
  }

  CouNumber x0 = info -> solution_ [ind],
            y0 = info -> solution_ [w -> Index ()];

  brpts = (double *) realloc (brpts, sizeof (double));

  // the best branching point for |x| is 0, as the two subproblems
  // will have exact convexifications
  *brpts = 0.;

  way = TWO_RAND; // don't really care which subtree to visit first

  // exact distance between current point and the two subsequent
  // convexifications
  return mymin (project (1., -1., 0., x0, y0, 0., COUENNE_INFINITY,  0, NULL, NULL),
		project (1., +1., 0., x0, y0, -COUENNE_INFINITY, 0., 0, NULL, NULL));
}
