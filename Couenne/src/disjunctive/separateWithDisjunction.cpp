/*
 * Name:    tightenWithDisjunction.cpp
 * Author:  Pietro Belotti
 * Purpose: preprocessing for disjunction -- return infeasible for
 *          none, one, or both sides of the disjunction
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"

/// generate one disjunctive cut from one CGLP
int CouenneDisjCuts::separateWithDisjunction (std::pair <OsiCuts *, OsiCuts *> &disj, 
					     OsiSolverInterface &si, 
					     OsiCuts &cs, 
					     const CglTreeInfo &info) const {

  return COUENNE_FEASIBLE;
}
