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
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg) {

  int w_ind = w         -> Index (),
      x_ind = argument_ -> Index ();

  expression *lbe, *ube;
  argument_ -> getBounds (lbe, ube);

  CouNumber l = (*lbe) (),
            u = (*ube) ();

  delete lbe;
  delete ube;

  bool cLeft  = !chg || (chg [x_ind].lower != UNCHANGED) || cg -> isFirst (),
       cRight = !chg || (chg [x_ind].upper != UNCHANGED) || cg -> isFirst ();

  // if l, u have the same sign, then w = x (l > 0) or w = -x (u < 0)

  if      (l >= -0) {if (cLeft)  cg -> createCut (cs, 0., 0, w_ind, 1., x_ind, -1.);}
  else if (u <=  0) {if (cRight) cg -> createCut (cs, 0., 0, w_ind, 1., x_ind, +1.);}
  else {

    // add two global cuts: w >= x and w >= -x
    if (cg -> isFirst ()) {
      cg -> createCut (cs, 0., +1, w_ind, 1., x_ind, -1.);
      cg -> createCut (cs, 0., +1, w_ind, 1., x_ind,  1.);
    }

    // otherwise check if at most one of the bounds is infinite: even
    // so, we can still add a plane, whose slope is 1 (x unbounded
    // from above) or -1 (from below)

    if (l > - COUENNE_INFINITY) {

      if (u < COUENNE_INFINITY) { // the upper approximation has slope other than -1, 1

	  CouNumber slope = (u+l) / (u-l);
	  // add an upper segment, which depends on the lower/upper bounds
	  if (cLeft || cRight) cg -> createCut (cs, -l*(slope+1.), -1, w_ind, 1., x_ind, -slope);
	}
	else // slope = 1
	  if (cLeft) cg -> createCut (cs, -2*l, -1, w_ind, 1., x_ind, -1.);
      }
      else if (u < COUENNE_INFINITY) // slope = -1
	if (cRight) cg -> createCut (cs, 2*u, -1, w_ind, 1., x_ind, 1.);
    }
}
