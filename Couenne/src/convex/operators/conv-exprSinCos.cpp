/*
 * Name:    conv-exprSinCos.C
 * Author:  Pietro Belotti
 * Purpose: convexification methods for sines and cosines
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>
#include <CouenneTypes.h>
#include <CouenneCutGenerator.h>
#include <exprSin.h>
#include <exprCos.h>
#include <exprAux.h>


// generate convexification cut for constraint w = sin (this)

void exprSin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  trigGenCuts (w, cs, cg, sin);
}


// generate convexification cut for constraint w = cos (this)

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  trigGenCuts (w, cs, cg, cos);
}
