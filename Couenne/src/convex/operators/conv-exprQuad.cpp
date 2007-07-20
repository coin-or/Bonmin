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
#include <exprBQuad.hpp>


/// Get lower and upper bound of an expression (if any)
void exprQuad::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLBQuad (this);
  ub = new exprUBQuad (this);

  //lb = new exprConst (-COUENNE_INFINITY);
  //ub = new exprConst  (COUENNE_INFINITY);
}


// generate equality between *this and *w
void exprQuad::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg, 
			     int wind, CouNumber lb, CouNumber ub) {

  // see if it is necessary to create/renew the alpha-convexification
  alphaConvexify (si);

  quadCuts (cs, cg);
}



