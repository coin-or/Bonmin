/*
 * Name:    conv-exprExp.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the exponential operator
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.hpp>
#include <exprExp.hpp>
#include <exprConst.hpp>
#include <exprAux.hpp>
#include <exprPow.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>


// generate convexification cut for constraint w = this

void exprExp::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs,  const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber l = (*le) (),
            u = (*ue) ();

  delete le;
  delete ue;

  int w_ind = aux       -> Index (),
      x_ind = argument_ -> Index ();

  bool cL = !chg || (chg [x_ind].lower != UNCHANGED) || cg -> isFirst (),
       cR = !chg || (chg [x_ind].upper != UNCHANGED) || cg -> isFirst ();

  // if bounds are very close, convexify with a single line

  if (fabs (u - l) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (u+l), ex0 = exp (x0);
    if (cL || cR) cg -> createCut (cs, ex0 * (1 - x0), 0, w_ind, 1., x_ind, - ex0);
    return;
  }

  CouNumber x = (cg -> isFirst ()) ? 
                 0 : powNewton ((*argument_) (), (*aux) (), exp, exp, exp);

  // upper segment

  if ((cL || cR) 
      && (u < log (COUENNE_INFINITY) ) 
      && (l > -    COUENNE_INFINITY / 1e4)) { // tame lower bound

    CouNumber expl     = exp (l),
              oppslope = (expl - exp (u)) / (u - l);

    cg -> createCut (cs, expl + oppslope*l, -1, 
		     w_ind, 1., 
		     x_ind, oppslope);
  }

  // add tangent points: first choose sampling points

  const CouNumber logMC = log (COU_MAX_COEFF);

  // add tangents with finite coefficients
  if (l < - logMC) l = - logMC;
  if (u >   logMC) u =   logMC;

  // approximate the exponential function from below
  cg -> addEnvelope (cs, +1, exp, exp, w_ind, x_ind, x, l, u, chg, true);
}
