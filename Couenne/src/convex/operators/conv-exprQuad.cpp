/*
 * Name:    conv-exprQuad.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: implementation of convexification methods for exprQuad
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>

#include <exprQuad.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

/// Get lower and upper bound of an expression (if any)
void exprQuad::getBounds (expression *&lb, expression *&ub) {

  expression *lbgrp, *ubgrp;

  /// compute lower/upper bound of nonlinear part
  exprGroup::getBounds (lbgrp, ubgrp);

  ///


}


// generate equality between *this and *w
void exprQuad::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg, 
			     int wind, CouNumber lb, CouNumber ub) {

  // see if it is necessary to create/renew the alpha-convexification
  alphaConvexify (si);

  // based on the information (dIndex_, dCoeffLo_, dCoeffUp_)
  // created/modified by alphaConvexify(), create convexification cuts
  // for this expression
  quadCuts (cs, cg);
}
