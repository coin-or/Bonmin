/*
 * Name:    exprMinMax.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of min and max operators
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>

#include <CouenneCutGenerator.hpp>
#include <CouenneTypes.h>
#include <exprMax.hpp>
#include <exprMin.hpp>

void exprMin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) 
{}


void exprMax::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) 
{}
