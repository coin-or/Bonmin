/*
 * Name:    branchExprPow.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for powers
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprPow.hpp>
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

  // TODO...

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
