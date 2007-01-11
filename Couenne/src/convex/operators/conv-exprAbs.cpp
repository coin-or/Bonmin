/*
 * Name:    conv-exprAbs.C
 * Author:  Pietro Belotti
 * Purpose: convexification methods for absolute value 
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>
#include <CouenneTypes.h>
#include <CouenneCutGenerator.h>
#include <exprAbs.h>
#include <exprAux.h>


// generate convexification cut for constraint w = |x|

void exprAbs::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  int w_ind = w         -> Index (),
      x_ind = argument_ -> Index ();

  // convexifying an expression w = |x| consists in adding two initial
  // constraints w >= x and w >= -x:

  if (cg -> isFirst ()) {

    cg -> addTangent (cs, w_ind, x_ind, 0, 0, CouNumber (-1.), +1);
    cg -> addTangent (cs, w_ind, x_ind, 0, 0, CouNumber  (1.), +1);
  }

  // add an upper segment, which depends on the lower/upper bounds

  expression *lbe, *ube;

  argument_ -> getBounds (lbe, ube);

  CouNumber l = (*lbe) (),
            u = (*ube) ();

  // if l and u have the same sign, then w = x (l > 0) or w = -x (u < 0)

  if      (l > 0) cg -> addTangent (cs, w_ind, x_ind, 0, 0, CouNumber (-1.), -1);
  else if (u < 0) cg -> addTangent (cs, w_ind, x_ind, 0, 0, CouNumber  (1.), -1);
  else {

    // otherwise check if one of the bounds is infinite, we can still
    // add a plane, whose slope will be one (unbounded from above) or
    // -1 (from below)

    if (l > - COUENNE_INFINITY + 1) {

      if (u < COUENNE_INFINITY - 1) // the upper approximation has slope other than -1, 1
	cg -> addSegment (cs, w_ind, x_ind, l, -l, u, u, -1);
      else // slope = 1
	cg -> addTangent (cs, w_ind, x_ind, l, -l, CouNumber (+1.), -1);
    }
    else if (u < COUENNE_INFINITY - 1) // slope = -1
      cg -> addTangent (cs, w_ind, x_ind, u, u, CouNumber (-1.), -1);
  }
}
