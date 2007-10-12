/*
 * Name:    CouenneSolverInterface.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Implementation of the OsiSolverInterface::resolve () method 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiClpSolverInterface.hpp>
#include <CouenneProblem.hpp>
#include <CouenneSolverInterface.hpp>

//#define DEBUG


/// Solve initial LP relaxation 
void CouenneSolverInterface::initialSolve () 
  {OsiClpSolverInterface::initialSolve ();}


/// defined in Couenne/src/convex/generateCuts.cpp
void sparse2dense (int, t_chg_bounds *, int *&, int &);


/// Resolve an LP relaxation after problem modification
void CouenneSolverInterface::resolve () {

  OsiClpSolverInterface::resolve ();
  return;

  int ncols = cutgen_ -> Problem () -> nVars ();

  cutgen_ -> Problem () -> update (getColSolution (),
				   getColLower    (),
				   getColUpper    ());

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  t_chg_bounds *chg_bds = new t_chg_bounds [ncols];

  OsiCuts cs;

  Bonmin::BabInfo *babInfo = dynamic_cast <Bonmin::BabInfo *> (getAuxiliaryInfo ());

  if (! (cutgen_ -> boundTightening (this, cs, chg_bds, babInfo))) {

#ifdef DEBUG
    printf ("#### BT says infeasible before re-solve\n");
#endif

    // how do we return infeasibility? TODO
  }

  int *changed = NULL, nchanged;
  sparse2dense (ncols, chg_bds, changed, nchanged);

  // change tightened bounds through OsiCuts
  if (nchanged)
    cutgen_ -> genColCuts (*this, cs, nchanged, changed);

  applyCuts (cs);

  // TODO: if NLP point available, add new cuts BEFORE resolving --
  // and decrease number of cutting plane iterations by one, to
  // balance it

  OsiClpSolverInterface::resolve ();
}
