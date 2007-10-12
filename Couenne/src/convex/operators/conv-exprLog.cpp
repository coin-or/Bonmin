/*
 * Name:    conv-exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the logarithm operator
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.hpp>
#include <exprLog.hpp>
#include <exprInv.hpp>
#include <exprPow.hpp>
#include <exprConst.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>


#define LOG_STEP 10
#define LOG_MININF 1e-50

// generate convexification cut for constraint w = this

void exprLog::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber l = (*le) (),
            u = (*ue) ();

  delete le;
  delete ue;

  int w_ind = aux       -> Index ();
  int x_ind = argument_ -> Index ();

  bool cL = !chg || (cg -> isFirst ()) || (chg [x_ind].lower != UNCHANGED),
       cR = !chg || (cg -> isFirst ()) || (chg [x_ind].upper != UNCHANGED);

  // if bounds are very close, convexify with a single line

  if ((fabs (u - l) < COUENNE_EPS) && (l > COUENNE_EPS)) {

    CouNumber x0 = 0.5 * (u+l);
    if (cL || cR) cg -> createCut (cs, log (x0) - 1, 0, 
				   w_ind, 1., x_ind, - 1/x0);
    return;
  }

  CouNumber x = (cg -> isFirst ()) ? 
                 1 : powNewton ((*argument_) (), (*aux) (), log, inv, oppInvSqr);

  // fix lower bound

  if (l < LOG_MININF) l = LOG_MININF;
  else   // lower segment (only put if lower bound is far enough from
	 // zero and upper is finite)
    if (u < COUENNE_INFINITY - 1) { 

      CouNumber dx   = u-l;
      CouNumber logu = log (u);
      CouNumber dw   = logu - log (l);

      if (cL || cR) 
	cg -> createCut (cs, dx*logu - u*dw, +1, w_ind, dx, x_ind, -dw);
    }

  // fix bound interval (unless you like infinite coefficients)

  if      (x < l) x = l;
  else if (x > u) x = u;

  if (u > 1e5 * log (COUENNE_INFINITY) - 1)
    u = x + (LOG_STEP << cg -> nSamples ());

  // add upper envelope

  cg -> addEnvelope (cs, -1, log, inv, w_ind, x_ind, x, l, u, chg, true);
}
