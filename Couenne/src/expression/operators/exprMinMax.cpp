/*
 * Name:    exprMinMax.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of min and max operators
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>

#include <CouenneCutGenerator.hpp>
#include <CouenneTypes.hpp>
#include <exprMax.hpp>
#include <exprMin.hpp>

void exprMin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) 
{}


void exprMax::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) 
{}
