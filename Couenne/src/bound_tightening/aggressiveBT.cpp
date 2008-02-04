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

bool CouenneProblem::aggressiveBT (t_chg_bounds *chg_bds, 
				   Bonmin::BabInfo * babInfo) const {
  int  ncols  = nVars ();
  bool retval = true;

  CouNumber
    *olb = new CouNumber [ncols],
    *oub = new CouNumber [ncols];

  // X is now the NLP solution, but in a low-dimensional space. We
  // have to get the corresponding point in higher dimensional space
  // through getAuxs()

  double *X = new double [ncols];
  CoinCopyN (babInfo -> nlpSolution (), nOrig_, X);
  getAuxs (X);

  // save current bounds
  CoinCopyN (Lb (), ncols, olb);
  CoinCopyN (Ub (), ncols, oub);

  // create a new, fictitious, bound bookkeeping structure
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    //    CouNumber cutoff = getCutOff ();
    int       objind = Obj (0) -> Body  () -> Index ();
    for (int i=0; i<nOrig_; i++)
      Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		      "   %2d %+20g [%+20g %+20g]\n",
		      i, X [i], Lb (i), Ub (i));
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		    "-------------\nAggressive BT. Current bound = %g, cutoff = %g, %d vars\n", 
		    Lb (objind), getCutOff (), ncols);
  }

  int improved, second, iter = 0;

  // Repeatedly fake tightening bounds on both sides of every variable
  // to concentrate around current NLP point.
  //
  // MAX_ABT_ITER is the maximum # of outer cycles. Each call to
  // fake_tighten in turn has an iterative algorithm for a
  // derivative-free, uni-dimensional optimization problem on a
  // monotone function.

  do {

    improved = 0;

    // scan all variables
    for (int i=0; i<ncols; i++) {

      int index = evalOrder (i);

      // AW: We only want to do the loop that temporarily changes
      // bounds around the NLP solution only for those points from the
      // NLP solution (no auxiliary vars)?

      // PBe: if we do want that, index should be initialized as i, as
      // evalOrder gives a variable index out of an array index.

      //      if (index < nOrig_) {

	// if (index == objind) continue; // don't do it on objective function

	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"------------- tighten left x%d\n", index);

	// tighten on left
	if ((X [index] >= Lb (index) + COUENNE_EPS)
	    && ((improved = fake_tighten (0, index, X, olb, oub, chg_bds, f_chg)) < 0)) {
	  retval = false;
	  break;
	}

	second = 0;

	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"------------- tighten right x%d\n", index);

	//Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"  ### tighten right\n");

	// tighten on right
	if ((X [index] <= Ub (index) - COUENNE_EPS)
	    && ((second = fake_tighten (1, index, X, olb, oub, chg_bds, f_chg) < 0))) {
	  retval = false;
	  break;
	}

	improved += second;
	//      }
    }
  } while (retval && improved && (iter++ < MAX_ABT_ITER));

  // store new valid bounds, or restore old ones if none changed
  CoinCopyN (olb, ncols, Lb ());
  CoinCopyN (oub, ncols, Ub ());

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"------------------\n");
    for (int i=0; i<ncols; i++)
      Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		      "   %2d %+20g %+20g  | %+20g\n", i, Lb (i), Ub (i), X [i]);
    if (!retval) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
				 "Couenne infeasible node from aggressive BT\n");
  }

  delete [] X;
  delete [] f_chg;
  delete [] olb;
  delete [] oub;

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    //    CouNumber cutoff = getCutOff ();
    int       objind = Obj (0) -> Body  () -> Index ();
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		    "-------------\ndone Aggressive BT. Current bound = %g, cutoff = %g, %d vars\n", 
		    Lb (objind), getCutOff (), ncols);
    for (int i=0; i<nOrig_; i++)
      Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		      "   %2d %+20g %+20g  | %+20g\n",
		      i, Lb (i), Ub (i), X [i]);
  }

  return retval;// && btCore (psi, cs, chg_bds, babInfo, true); // !!!
  //return retval && btCore (psi, cs, chg_bds, babInfo, true);
}
