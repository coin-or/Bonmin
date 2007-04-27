/*
 * Name:    conv-exprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the exponential operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprAux.h>
#include <exprPow.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// generate convexification cut for constraint w = this

void exprExp::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs,  const CouenneCutGenerator *cg) {
  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber x = (cg -> isFirst ()) ? 
                 0 : powNewton ((*argument_) (), (*aux) (), exp, exp, exp),
            l = (*le) (),
            u = (*ue) ();

  int w_ind = aux       -> Index (),
      x_ind = argument_ -> Index ();

  // upper segment

  if ((   u < log (COUENNE_INFINITY) - 1) 
      && (l > -    COUENNE_INFINITY  + 1)) {

    CouNumber expl     = exp (l),
              oppslope = (expl - exp (u)) / (u - l);

    cg -> createCut (cs, expl + oppslope*l, -1, 
		     w_ind, 1., 
		     x_ind, oppslope);
  }

  // add tangent points: first choose sampling points

  int ns = cg -> nSamples ();

  // add tangents with finite coefficients
  if (l < log (COUENNE_EPS))      l = x - ns;
  if (u > log (COUENNE_INFINITY)) u = x + ns;

  // approximate the exponential function from below
  cg -> addEnvelope (cs, +1, exp, exp, w_ind, x_ind, x, l, u, true);

  delete le;
  delete ue;
}
