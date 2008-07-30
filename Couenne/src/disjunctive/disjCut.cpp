/*
 * Name:    disjCut.cpp
 * Author:  Pietro Belotti
 * Purpose: generate one disjunctive cut based on a single disjunction
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"

/// generate one disjunctive cut from one CGLP
int CouenneDisjCuts::generateDisjCut (std::pair <OsiCuts *, OsiCuts *> &disj, 
				      OsiSolverInterface &si, 
				      OsiCuts &cs, 
				      const CglTreeInfo &info) const {
  OsiCuts
    *left  = disj.first,
    *right = disj.second;

  // merge si+left and si+right in CGLP

  return COUENNE_FEASIBLE;
}
