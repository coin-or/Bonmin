/*
 * Name:    conv-exprAbs.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for |f(x)|
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

  // convexifying an expression w = |x| consists in adding two
  // global cuts: w >= x and w >= -x

  if (cg -> isFirst ()) {

    cg -> createCut (cs, 0., +1, w_ind, 1., x_ind, -1.);
    cg -> createCut (cs, 0., +1, w_ind, 1., x_ind,  1.);
  }

  // add an upper segment, which depends on the lower/upper bounds

  expression *lbe, *ube;

  argument_ -> getBounds (lbe, ube);

  CouNumber l = (*lbe) (),
            u = (*ube) ();

  // if l, u have the same sign, then w = x (l > 0) or w = -x (u < 0)

  if      (l >= 0) cg -> createCut (cs, 0., -1, w_ind, 1., x_ind, -1.);
  else if (u <= 0) cg -> createCut (cs, 0., -1, w_ind, 1., x_ind, +1.);
  else {

      // otherwise check if only one of the bounds is infinite: if so,
      // we can still add a plane, whose slope will be one (if x is
      // unbounded from above) or -1 (from below)

      if (l > - COUENNE_INFINITY + 1) {

	if (u < COUENNE_INFINITY - 1) { // the upper approximation has slope other than -1, 1

	  CouNumber slope = (u+l) / (u-l);
	  cg -> createCut (cs, -l*(slope+1.), -1, w_ind, 1., x_ind, -slope);
	}
	else // slope = 1
	  cg -> createCut (cs, -2*l, -1, w_ind, 1., x_ind, -1.);
      }
      else if (u < COUENNE_INFINITY - 1) // slope = -1
	cg -> createCut (cs, 2*u, -1, w_ind, 1., x_ind, 1.);
    }

  delete lbe;
  delete ube;
}
