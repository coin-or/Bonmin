/*
 * Name:    conv-exprQuad.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: implementation of convexification methods for exprQuad
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"

#include "exprAux.hpp"
#include "exprQuad.hpp"
#include "exprBQuad.hpp"
#include "CouenneCutGenerator.hpp"

/// Get lower and upper bound of an expression (if any)
void exprQuad::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLBQuad (this);
  ub = new exprUBQuad (this);
}


// generate equality between *this and *w
void exprQuad::generateCuts (expression *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg, 
			     int wind, CouNumber lb, CouNumber ub) {

  if ((!(cg -> isFirst ())) &&                    // unless a convexification was never created,
      (fabs ((*this) () - (*w) ()) < COUENNE_EPS) // do we really need a convexification cut?
      || !alphaConvexify (cg -> Problem (), si))  // ... or a new alpha-convexification?
    return;

  // generate linear cuts for convex quadratic [upper|lower]-envelope
  // of this expression
  quadCuts (w, cs, cg);
}
