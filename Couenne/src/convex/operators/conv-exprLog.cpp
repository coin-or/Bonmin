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
#define LOG_SCALE 1e-20

// generate convexification cut for constraint w = this

void exprLog::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {
  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber// w = (*aux) (),
            x = (*argument_) (),
            l = (*le) (),
            u = (*ue) ();

  OsiRowCut *cut;

  int w_ind = aux       -> Index ();
  int x_ind = argument_ -> Index ();

  // fix lower bound

  if (l < COUENNE_EPS) 
    l = COUENNE_EPS;
  else   // lower segment (only put if lower bound is far enough from
	 // zero and upper is finite
    if (u < COUENNE_INFINITY - 1) { 

      CouNumber dx   = u-l;
      CouNumber logu = log (u);
      CouNumber dw   = logu - log (l);

      if ((cut = cg -> createCut (dx*logu - u*dw, +1, w_ind, dx, 
				  x_ind, -dw)))
	cs.insert (cut);
    }

  // fix bound interval (unless you like infinite coefficients)

  int ns = cg -> nSamples ();

  if (u > COUENNE_INFINITY - 1)
    u = l + (LOG_STEP << ns);

  if (l < 2 * COUENNE_EPS)
    l = LOG_SCALE;

  if (x <= COUENNE_EPS)
    l = LOG_SCALE;

  // add upper envelope

  cg -> addEnvelope (cs, -1, log, inv, w_ind, x_ind, x, l, u, true);

  delete le;
  delete ue;
}
