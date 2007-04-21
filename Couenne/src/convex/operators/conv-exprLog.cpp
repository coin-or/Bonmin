/*
 * Name:    conv-exprLog.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the logarithm operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprLog.h>
#include <exprConst.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


#define LOG_STEP 10
#define LOG_MININF 1e-300

// generate convexification cut for constraint w = this

void exprLog::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber x = (*argument_) (),
            l = (*le) (),
            u = (*ue) ();

  int w_ind = aux       -> Index ();
  int x_ind = argument_ -> Index ();

  // fix lower bound

  if (l < LOG_MININF) l = LOG_MININF;
  else   // lower segment (only put if lower bound is far enough from
	 // zero and upper is finite)
    if (u < COUENNE_INFINITY - 1) { 

      CouNumber dx   = u-l;
      CouNumber logu = log (u);
      CouNumber dw   = logu - log (l);

      cg -> createCut (cs, dx*logu - u*dw, +1, w_ind, dx, x_ind, -dw);
    }

  // fix bound interval (unless you like infinite coefficients)

  if      (x < l) x = l;
  else if (x > u) x = u;

  if (u > COUENNE_INFINITY - 1)
    u = x + (LOG_STEP << cg -> nSamples ());

  // add upper envelope

  cg -> addEnvelope (cs, -1, log, inv, w_ind, x_ind, x, l, u, true);

  delete le;
  delete ue;
}
