/*
 * Name:    branchExprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for logarithms
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprLog.h>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprLog::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {
  ind = -1;
  return 0.;

  /*  ind = argument_ -> Index ();

  if (ind < 0) {
    printf ("argument of exprAbs has negative index\n");
    exit (-1);
  }

  brpts = (double *) malloc (sizeof (double));
  *brpts = 0.;
  way = TWO_LEFT;

  return mymin (project (1., -1., 0., x0, y0, 0., COUENNE_INFINITY,  0, NULL, NULL),
  project (1., +1., 0., x0, y0, -COUENNE_INFINITY, 0., 0, NULL, NULL));*/
}

/*
/// distance covered by current point if branching rule applied to this expression
double exprLog::BranchGain (expression *w, const OsiBranchingInformation *info) {

  int xi = argument_ -> Index (),
      wi = w         -> Index ();

  if (xi < 0 || wi < 0) {
    printf ("negative indices inbranchExprAbs\n");
    exit (-1);
  }

  CouNumber x0 = expression::Variable (xi),
            y0 = expression::Variable (wi);

  return 0;
  //  return mymin (project (1, -1, 0, x0, y0, 0, COUENNE_INFINITY,  0, NULL, NULL),
  //		project (1, +1, 0, x0, y0, -COUENNE_INFINITY, 0, 0, NULL, NULL));
}

/// branching object best suited for this expression
OsiBranchingObject *exprLog::BranchObject (expression *w, const OsiBranchingInformation *) {

  // branching once on 0 is sufficient for expressions w=|x| 
  //  return new CouenneBranchingObject (Argument (), 0);
  return NULL;
}
*/
