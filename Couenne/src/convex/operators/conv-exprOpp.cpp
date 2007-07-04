/*
 * Name:    conv-exprOpp.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to convexify opposite of expressions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOpp.hpp>
#include <exprConst.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>


// generate convexification cut for constraint w = - x

void exprOpp::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) {
  // easy... 

  if (cg -> isFirst ())
    cg -> createCut (cs, 0., 0, w->Index (), 1., argument_->Index (), 1.);
}
