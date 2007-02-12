/*
 * Name: generateCuts.cpp
 * Author: Pietro Belotti
 * Purpose: the generateCuts() method of the convexification class CouenneCutGenerator
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>


// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si, 
					OsiCuts &cs, 
					const CglTreeInfo info) const {
  int ncols = si.getNumCols ();

  CouNumber *x, *l, *u;

  if (firstcall_) {

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to 

    // For each auxiliary variable replacing the original constraints,
    // check if corresponding bounds are violated, and add cut to cs

    int nnlc = problem_ -> nNLCons ();

    for (int i=0; i<nnlc; i++) {

      CouenneConstraint *con = problem_ -> NLCon (i);

      // for constraint lb <= w <= ub, compute actual values of lb, ub

      CouNumber lb = con -> Lb () -> Value ();
      CouNumber ub = con -> Ub () -> Value ();

      // if there exists violation, add constraint

      OsiRowCut *orc   = new OsiRowCut;
      CouNumber *coeff = new CouNumber [1];
      int       *index = new int       [1];

      coeff [0] = 1;
      index [0] = con -> Body () -> Index ();

      orc -> setRow (1, index, coeff);

      if (lb > - COUENNE_INFINITY + 1) orc -> setLb (lb);
      if (ub <   COUENNE_INFINITY - 1) orc -> setUb (ub);

      cs.insert (orc);
    }
  }
  else {

    // Retrieve, from si, value and bounds of all variables, if not
    // firstcall, otherwise only those of the original ones Update
    // expression structure with x, l, u

    const CouNumber *xc = si.getColSolution (), 
      *lc = si.getColLower (),
      *uc = si.getColUpper ();

    x = new CouNumber [ncols];
    l = new CouNumber [ncols];
    u = new CouNumber [ncols];

    for (int i=ncols; i--;) {

      x [i] = xc [i];
      l [i] = lc [i];
      u [i] = uc [i];
    }

    problem_ ->  update (x,l,u);
    //    expression:: update (x,l,u);
  }

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  for (int i=0; i<problem_ -> nAuxs (); i++)
    problem_ -> Aux (i) -> generateCuts (si, cs, this);

  // end of generateCuts

  if (cs.sizeRowCuts ())
    printf ("Couenne: %d convexifier cuts\n", cs.sizeRowCuts ());

  if (firstcall_) 
    firstcall_ = false;
  else {
    delete [] x;
    delete [] l;
    delete [] u;
  }
}
