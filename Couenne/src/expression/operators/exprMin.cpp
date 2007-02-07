/*
 * Name:    exprMin.C
 * Author:  Pietro Belotti
 * Purpose: convexification methods for min operator 
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>
#include <CouenneCutGenerator.h>
#include <exprMin.h>

void exprMin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) 
{}
