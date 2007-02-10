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

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// generate convexification cut for constraint w = this

void exprExp::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs,  const CouenneCutGenerator *cg) {

  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber x = (*argument_) (),
            l = (*le) (),
            u = (*ue) ();

  int w_ind = aux       -> Index (),
      x_ind = argument_ -> Index ();

  OsiRowCut *cut;

  int j = cs.sizeRowCuts ();

  // upper segment

  if ((   u < log (COUENNE_INFINITY) - 1) 
      && (l > -    COUENNE_INFINITY  + 1)) {

    CouNumber expl     = exp (l),
              oppslope = (expl - exp (u)) / (u - l);

    if ((cut = cg -> createCut (expl + oppslope*l, -1, 
				w_ind, CouNumber (1.), 
				x_ind, oppslope)))
      cs.insert (cut);
  }

  // lower convexification: start with trivial envelope w >= 0

  if ((cut = cg -> createCut (CouNumber (0.), +1, w_ind, CouNumber (1.))))
    cs.insert (cut);

  // add tangent points: first choose sampling points

  int ns = cg -> nSamples ();

  // fix bounds to get finite coefficients

  CouNumber fact = 2 * ns;

  if (x > COUENNE_EPS) {
    l = x / fact;
    u = x * fact;
  }
  else if (x < -COUENNE_EPS) {
    l = x * fact;
    u = x / fact;
  }
  else {
    l = - ns;
    u = ns;
  }
    
  // approximate the exponential function from below
  cg -> addEnvelope (cs, +1, exp, exp, w_ind, x_ind, x, l, u, true);

  delete le;
  delete ue;
}
