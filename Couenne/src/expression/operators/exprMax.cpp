/*
 * Name:    exprMax.cpp
 * Author:  Pietro Belotti
 * Purpose: definition of max
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprMax.h>
#include <exprCopy.h>

// Get lower and upper bound of an expression (if any)
/*
void exprMax::getBounds (expression *&lb, expression *&ub) {

  expression **all = new expression * [nargs_];
  expression **alu = new expression * [nargs_];

  for (register int i=0; i<nargs_; i++) {

    all [i] = new exprCopy (arglist_ [i]); // i is even
    alu [i] = new exprCopy (arglist_ [i]);

    i++; // now i is odd

    arglist_ [i] -> getBounds (all [i], alu [i]);
  }

  lb = new exprMax (all, nargs_);
  ub = new exprMax (alu, nargs_);
}
*/
// generate equality between *this and *w
void exprMax::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) 
{}
