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

#define MAX_ABT_ITER 8  // max # aggressive BT iterations


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

  const double *X = babInfo->nlpSolution(); // length nOrig()

  // save current bounds
  CoinCopyN (lb, ncols, olb);
  CoinCopyN (ub, ncols, oub);

  // create new, fictitious, bounds
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    CouNumber cutoff = problem_ -> getCutOff ();
    int       objind = problem_ -> Obj (0) -> Body  () -> Index ();
    for (int i=0; i<problem_ -> nOrig(); i++)
      Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		      "   %2d %+20g %+20g  | %+20g\n",
		      i, lb [i], ub [i], X [i]);
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		    "-------------\nAggressive BT. Current bound = %g, cutoff = %g, %d vars\n", 
		    lb [objind], cutoff, ncols);
  }

  int improved = 0, second, iter = 0;

  // Repeatedly fake tightening bounds on both sides of every variable
  // to concentrate around current NLP point.
  //
  // MAX_ABT_ITER is the maximum # of outer cycles. Each call to
  // fake_tighten in turn has an iterative algorithm for a
  // derivative-free, uni-dimensional optimization problem on a
  // monotone function.

  do {

    improved = 0;  // Pietro: This statement was inside for loop - correct here?

    // scan all variables
    for (int i=0; i<ncols; i++) {

      int index = problem_ -> evalOrder (i);

      // AW: We only want to do the loop that temporarily changes bounds
      // around the NLP solution only for those points from the NLP
      // solution (no auxiliary vars)?
      if (index < problem_ -> nOrig()) {

	// if (index == objind) continue; // don't do it on objective function

	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"x_%03d:-----------------------------\n  ### tighten left\n", index);

	// tighten on left
	if ((X [index] >= lb [index] + COUENNE_EPS)
	    && ((improved = fake_tighten (psi, cs, babInfo, 
					  0, index, X, olb, oub, chg_bds, f_chg)) < 0)) {
	  retval = false;
	  break;
	}

	second = 0;

	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"  ### tighten right\n");

	// tighten on right
	if ((X [index] <= ub [index] - COUENNE_EPS)
	    && ((second = fake_tighten (psi, cs, babInfo, 
				      1, index, X, olb, oub, chg_bds, f_chg) < 0))) {
	  retval = false;
	  break;
	}

	improved += second;
      }
    }

  } while (retval && improved && (iter++ < MAX_ABT_ITER));

  // store new valid bounds into problem, or restore old ones if none changed
  CoinCopyN (olb, ncols, lb);
  CoinCopyN (oub, ncols, ub);

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"------------------\n");
    for (int i=0; i<ncols; i++)
      Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"   %2d %+20g %+20g  | %+20g\n", i, lb [i], ub [i], X [i]);
    if (!retval) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"#### infeasible node from aggressive BT\n");
  }

  delete [] f_chg;
  delete [] olb;
  delete [] oub;

  return retval;// && btCore (psi, cs, chg_bds, babInfo, true); // !!!
  //return retval && btCore (psi, cs, chg_bds, babInfo, true);
}
