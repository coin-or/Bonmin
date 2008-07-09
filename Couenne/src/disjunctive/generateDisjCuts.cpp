/*
 * Name:    generateDisjCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: separation method for disjunctive cuts
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinTime.hpp"

#include "CouenneDisjCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

/// generate disjunctive cuts
void CouenneDisjCuts::generateCuts (const OsiSolverInterface &si, 
				    OsiCuts &cs, 
				    const CglTreeInfo info) const {

  double time  = CoinCpuTime ();
  int    ncuts = 0;

  // consider problem Ax <= a_0, x in [l,u]
  //
  // get set of disjunctions (in the form of OsiBranchingObjects) in
  // number limited by depth, problem size, and sorted according to
  // branching method (HotInfo if strong branching?)
  //
  // for each disjunction (x_i <= or >= x_i^d) 
  //
  //   1a) apply left  disj., get cuts Bx <= b_0, x in [l_1,u_1]
  //   1b) apply right disj., get cuts Cx <= c_0, x in [l_2,u_2]
  //
  //   2a) if both disjunctions above are infeasible, return
  //       infeasible
  //
  //   2b) if one disjunction is infeasible, apply bound tightening
  //       and restart
  //
  //   2c) generate CGLP:
  //         consider subset of Ax <= a_0: 
  //           a) active only, or
  //           b) {root linearization cuts} union {currently active}, or
  //           c) all cuts (!)
  //
  //   3) solve CGLP
  //
  //   4) add corresponding cut

  if (info.level < 0) 
    nrootcuts_ = ncuts;
  ntotalcuts_ += ncuts;

  septime_ += (CoinCpuTime () - time);
}
