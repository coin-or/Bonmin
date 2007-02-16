/*
 * Name:    exprMinMax.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of min and max operators
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>
#include <CouenneCutGenerator.h>
#include <exprMax.h>
#include <exprMin.h>

void exprMin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) 
{}


void exprMax::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) 
{}
