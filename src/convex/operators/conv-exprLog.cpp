/*
 * Name:    conv-exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the logarithm operator
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"
#include "exprLog.hpp"
#include "exprPow.hpp"
#include "exprAux.hpp"

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"


#define LOG_STEP 10
#define LOG_MININF 1e-50

// generate convexification cut for constraint w = this

void exprLog::generateCuts (expression *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  int 
    w_ind = aux       -> Index (), //   dependent variable's index
    x_ind = argument_ -> Index (); // independent variable's index

  bool changed = 
    !chg || (cg -> isFirst ()) || 
    (chg [x_ind].lower() != t_chg_bounds::UNCHANGED) ||
    (chg [x_ind].upper() != t_chg_bounds::UNCHANGED);

  CouNumber l, u;
  argument_ -> getBounds (l, u);

  // if bounds are very close, convexify with a single line
  if ((fabs (u - l) < COUENNE_EPS) && (l > COUENNE_EPS)) {

    CouNumber x0 = 0.5 * (u+l);
    if (changed) cg -> createCut (cs, log (x0) - 1, 0, w_ind, 1., x_ind, - 1/x0);
    return;
  }

  // fix lower bound
  if (l < LOG_MININF) l = LOG_MININF;
  else   // lower segment (if l is far from zero and u finite)
    if ((u < COUENNE_INFINITY) && changed) { 

      CouNumber dx   = u-l;
      CouNumber logu = log (u);
      CouNumber dw   = logu - log (l);

      cg -> createCut (cs, dx*logu - u*dw, +1, w_ind, dx, x_ind, -dw);
    }

  // pick tangent point closest to current point (x0,y0)
  CouNumber x = (cg -> isFirst ()) ? 
                 1 : powNewton ((*argument_) (), (*aux) (), log, inv, oppInvSqr);

  // check if outside interval
  if      (x < l) x = l;
  else if (x > u) x = u;

  // fix bound interval (unless you like infinite coefficients)
  if (u > 1e5 * log (COUENNE_INFINITY) - 1)
    u = x + (LOG_STEP << cg -> nSamples ());

  // add upper envelope
  cg -> addEnvelope (cs, -1, log, inv, w_ind, x_ind, x, l, u, chg, true);
}
