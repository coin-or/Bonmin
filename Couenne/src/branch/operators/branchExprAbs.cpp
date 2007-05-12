/*
 * Name:    branchExprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for exprAbs
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprAbs.h>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneBranchingObject.hpp>
#include <CouenneObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprAbs::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {
  ind = argument_ -> Index ();

  CouNumber x0 = info -> solution_ [ind],
            y0 = info -> solution_ [w -> Index ()];

  if (ind < 0) {
    printf ("argument of exprAbs has negative index\n");
    exit (-1);
  }

  brpts = (double *) realloc (brpts, sizeof (double));
  *brpts = 0.;
  way = TWO_RAND; // don't really care which subtree to visit first

  return mymin (project (1., -1., 0., x0, y0, 0., COUENNE_INFINITY,  0, NULL, NULL),
		project (1., +1., 0., x0, y0, -COUENNE_INFINITY, 0., 0, NULL, NULL));
}

/*
/// distance covered by current point if branching rule applied to this expression
double exprAbs::BranchGain (expression *w, const OsiBranchingInformation *info) {

  int xi = argument_ -> Index (),
      wi = w         -> Index ();

  if (xi < 0 || wi < 0) {
    printf ("negative indices inbranchExprAbs\n");
    exit (-1);
  }

  CouNumber x0 = expression::Variable (xi),
            y0 = expression::Variable (wi);

  return
}

/// branching object best suited for this expression
OsiBranchingObject *exprAbs::BranchObject (expression *w, const OsiBranchingInformation *) {

  // branching once on 0 is sufficient for expressions w=|x| 
  return new CouenneBranchingObject (Argument (), 0);
}
*/
