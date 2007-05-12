/*
 * Name:    branchExprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for powers
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprPow.h>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>
#include <projections.h>


/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprPow::selectBranch (expression *w, 
				 const OsiBranchingInformation *info,
				 int &ind, 
				 double * &brpts, 
				 int &way) {
  ind = -1;
  return 0.;

  //  ind = argument_ -> Index ();

  /*  if (ind < 0) {
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
double exprPow::BranchGain (expression *w, const OsiBranchingInformation *info) {

  int xi = arglist_ [0] -> Index (),
      wi = w            -> Index ();

  if (xi < 0 || wi < 0) {
    printf ("negative indices in branchExprPow\n");
    exit (-1);
  }

  if (arglist_ [1] -> Type () != CONST) {
    printf ("non-constant exponent in branchExprPow\n");
    exit (-1);
  }

  CouNumber x0 = expression::Variable (xi),
            y0 = expression::Variable (wi),
            k  = arglist_ [1] -> Value ();

  return 0.;

  //  return mymin (project (1, -1, 0, x0, y0, 0, COUENNE_INFINITY,  0, NULL, NULL),
  //		project (1, +1, 0, x0, y0, -COUENNE_INFINITY, 0, 0, NULL, NULL));
}

/// branching object best suited for this expression
OsiBranchingObject *exprPow::BranchObject (expression *w, const OsiBranchingInformation *) {

  // branching once on 0 is sufficient for expressions w=|x| 
  //  return new CouenneBranchingObject (Argument (), 0);
  return NULL;
}
*/
