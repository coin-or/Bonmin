/*
 * Name:    aggressiveBT.cpp
 * Author:  Pietro Belotti
 * Purpose: aggressive bound tightening -- fake bounds in variables to
 *          exclude parts of the solution space through fathoming on
 *          bounds/infeasibility
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "BonBabInfos.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

//#define DEBUG

#define MAX_ABT_ITER 8  // max # aggressive BT iterations

// core of the bound tightening procedure
bool btCore (const CouenneCutGenerator *cg,
	     const OsiSolverInterface *psi,
	     OsiCuts &cs, 
	     t_chg_bounds *chg_bds, 
	     Bonmin::BabInfo * babInfo,
	     bool serious);


// single fake tightening. Return
//
// -1   if infeasible
//  0   if no improvement
// +1   if improved
int fake_tighten (const CouenneCutGenerator *cg,
		  const OsiSolverInterface *psi,
		  OsiCuts &cs,
		  Bonmin::BabInfo * babInfo,

		  char direction,  // 0: left, 1: right
		  int index,       // index of the variable tested
		  CouNumber *olb,  // cur. lower bound
		  CouNumber *oub,  //      upper
		  t_chg_bounds *chg_bds,
		  t_chg_bounds *f_chg);


// Aggressive Bound Tightening: for each variable, fake new bounds
// [l,b] or [b,u] and apply bound tightening. If the interval is
// fathomed on bounds or on infeasibility, the complementary bound
// interval is a valid tightening.

bool CouenneCutGenerator::aggressiveBT (const OsiSolverInterface *psi,
					OsiCuts &cs, 
					t_chg_bounds *chg_bds, 
					Bonmin::BabInfo * babInfo) const {
  int  ncols  = problem_ -> nVars ();
  bool retval = true;

  CouNumber
    *olb = new CouNumber [ncols],
    *oub = new CouNumber [ncols],
    *lb  = problem_ -> Lb (),
    *ub  = problem_ -> Ub ();

  const double *X = psi -> getColSolution ();

  // save current bounds
  CoinCopyN (lb, ncols, olb);
  CoinCopyN (ub, ncols, oub);

  // create new, fictitious, bounds
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

#ifdef DEBUG
  CouNumber cutoff = problem_ -> getCutOff ();
  int       objind = problem_ -> Obj (0) -> Body  () -> Index ();
  for (int i=0; i<ncols; i++)
    printf ("   %2d %+20g %+20g  | %+20g\n", i, lb [i], ub [i], X [i]);
  printf ("-------------\nAggressive BT. Current bound = %g, cutoff = %g, %d vars\n", 
	  lb [objind], cutoff, ncols);
#endif

  int improved = 0, second, iter = 0;

  // Repeatedly fake tightening bounds on both sides of every variable
  // to concentrate around current NLP point.
  //
  // MAX_ABT_ITER is the maximum # of outer cycles. Each call to
  // fake_tighten in turn has an iterative algorithm for a
  // derivative-free, uni-dimensional optimization problem on a
  // monotone function.

  do {

    // scan all variables
    for (int i=0; i<ncols; i++) {

      int index = problem_ -> evalOrder (i);

      // if (index == objind) continue; // don't do it on objective function

      improved = 0;

#ifdef DEBUG
      printf ("x_%03d:-----------------------------\n  ### tighten left\n", index);
#endif

      // tighten on left
      if ((X [index] >= lb [index] + COUENNE_EPS)
	  && ((improved = fake_tighten (this, psi, cs, babInfo, 
					0, index, olb, oub, chg_bds, f_chg)) < 0)) {
	retval = false;
	break;
      }

#ifdef DEBUG
      printf ("  ### tighten right\n");
#endif

      // tighten on right
      if ((X [index] <= ub [index] - COUENNE_EPS)
	  && ((second = fake_tighten (this, psi, cs, babInfo, 
				      1, index, olb, oub, chg_bds, f_chg) < 0))) {
	retval = false;
	break;
      }

      improved += second;
    }

  } while (retval && improved && (iter < MAX_ABT_ITER));

  // store new valid bounds into problem, or restore old ones if none changed
  CoinCopyN (olb, ncols, lb);
  CoinCopyN (oub, ncols, ub);

#ifdef DEBUG
  printf ("------------------\n");
  for (int i=0; i<ncols; i++)
    printf ("   %2d %+20g %+20g  | %+20g\n", i, lb [i], ub [i], X [i]);
  if (!retval) printf ("#### infeasible node from aggressive BT\n");
#endif

  delete [] f_chg;
  return retval;// && btCore (this, psi, cs, chg_bds, babInfo, true); // !!!
  //return retval && btCore (this, psi, cs, chg_bds, babInfo, true);
}
