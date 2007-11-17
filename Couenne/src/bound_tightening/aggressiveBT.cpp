/*
 * Name:    aggressiveBT.cpp
 * Author:  Pietro Belotti
 * Purpose: aggressive bound tightening -- fake bounds in variables to
 *          exclude parts of the solution space through fathoming on
 *          bounds/infeasibility
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "BonBabInfos.hpp"

// max # bound tightening iterations
#define MAX_BT_ITER 8
#define THRES_IMPROVED 0

//#define DEBUG

// the larger, the more conservative. Must be > 0
#define AGGR_MUL 2

// the smaller, the more conservative. Must be > 1
#define AGGR_DIV 20

//#define DEBUG

// core of the bound tightening procedure

bool btCore (const CouenneCutGenerator *cg,
	     const OsiSolverInterface *psi,
	     OsiCuts &cs, 
	     t_chg_bounds *chg_bds, 
	     Bonmin::BabInfo * babInfo,
	     bool serious);


void fictBounds (CouNumber  x,
		 CouNumber  lb,   CouNumber  ub, 
		 CouNumber &left, CouNumber &right) {

  if (lb < -COUENNE_INFINITY) {
    if (ub >  COUENNE_INFINITY) { // ]-inf,+inf[

      if (fabs (x) < COUENNE_EPS) left = - (right = AGGR_MUL);
      else                        left = - (right = AGGR_MUL * fabs (x));

    } else { // ]-inf,u]

      if      (x < -COUENNE_EPS) {left = AGGR_MUL * x;  right = CoinMin (0., (x+ub)/2);}
      else if (x >  COUENNE_EPS) {left = 0;             right = x + (ub-x)/AGGR_DIV;}
      else                       {left = -AGGR_MUL;     right =         ub/AGGR_DIV;}
    }
  }
  else {
    if (ub >  COUENNE_INFINITY) { // [l,+inf[

      if      (x < -COUENNE_EPS) {left = x - (x-lb)/AGGR_DIV;   right = 0;}
      else if (x >  COUENNE_EPS) {left = CoinMax (0.,(x+lb)/2); right = AGGR_MUL * x;}
      else                       {left = lb/AGGR_DIV;           right = AGGR_MUL;}
    } else { // [l,u]

      left  = x - (x-lb)/AGGR_DIV;
      right = x + (ub-x)/AGGR_DIV;
    }
  }
}


// Aggressive Bound Tightening
//
// For each variable, fake new bounds [l,b] or [b,u] and apply bound
// tightening. If the interval is fathomed on bounds or on
// infeasibility, the complementary bound interval is a valid
// tightening.

bool CouenneCutGenerator::aggressiveBT (const OsiSolverInterface *psi,
					OsiCuts &cs, 
					t_chg_bounds *chg_bds, 
					Bonmin::BabInfo * babInfo) const {
  int ncols = psi -> getNumCols ();

  bool retval = false;

  // save current bounds

  CouNumber
    *olb = new CouNumber [ncols],
    *oub = new CouNumber [ncols],
    *lb  = problem_ -> Lb (),
    *ub  = problem_ -> Ub (),
    cutoff = problem_ -> getCutOff ();

  const double *X = psi -> getColSolution ();

  CoinCopyN (lb, ncols, olb);
  CoinCopyN (ub, ncols, oub);

  int 
    objsense = problem_ -> Obj (0) -> Sense (),
    objind   = problem_ -> Obj (0) -> Body  () -> Index ();

  assert (objind >= 0);

  // create new bounds

  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

#ifdef DEBUG
  for (int i=0; i<ncols; i++)
    printf ("   %2d %+20g %+20g  | %+20g\n", i, lb [i], ub [i], X [i]);
  printf ("-------------\nAggressive BT. Current bound = %g, cutoff = %g\n", lb [objind], cutoff);
#endif

  for (int i=0; i<ncols; i++) 

    if (i != objind) { // don't do it on objective function

      CouNumber flb, fub;
      bool 
	feasLeft  = true, 
	feasRight = true,
	betterbds = false;

      fictBounds (X [i], olb [i], oub [i], flb, fub);

#ifdef DEBUG
      printf ("fake bounds %d: [%g | %g %g | %g], x = %g\n", 
	      i, olb [i], flb, fub, oub [i], X [i]);
#endif

      if (flb > lb [i] + COUENNE_EPS) {

	// fake left interval
	ub [i] = flb;
	feasLeft = btCore (this, psi, cs, f_chg, babInfo, false);
	betterbds = (objsense == MINIMIZE) && (lb [objind] > cutoff) || 
                    (objsense == MAXIMIZE) && (ub [objind] < cutoff);
	CoinCopyN (chg_bds, ncols, f_chg);
	CoinCopyN (olb, ncols, lb);
	CoinCopyN (oub, ncols, ub);

	if (!feasLeft || betterbds) { // left interval is feasible, check 

#ifdef DEBUG
	  printf ("[%2d] prune left  [%g --> %g, %g] <%d> -- %g %g\n", 
		  i, olb [i], flb, oub [i], feasLeft, lb [objind], cutoff);
#endif

	  // left interval can be left out
	  feasLeft = false;
	  olb [i] = lb [i] = flb;
	  chg_bds [i].setLower (t_chg_bounds::CHANGED);
	  if (!(btCore (this, psi, cs, chg_bds, babInfo, true)))
	    goto end_procedure;
	}
      }

      if (fub < ub [i] - COUENNE_EPS) {

	// fake right interval
	lb [i] = fub;
	feasRight = btCore (this, psi, cs, f_chg, babInfo, false);
	betterbds = (objsense == MINIMIZE) && (lb [objind] > cutoff) || 
                    (objsense == MAXIMIZE) && (ub [objind] < cutoff);
	CoinCopyN (chg_bds, ncols, f_chg);
	CoinCopyN (olb, ncols, lb);
	CoinCopyN (oub, ncols, ub);

	if (!feasRight || betterbds) { // left interval is feasible, check 

#ifdef DEBUG
	  printf ("[%2d] prune right [%g, %g <-- %g] <%d> -- %g %g\n", 
		  i, olb [i], fub, oub [i], feasRight, lb [objind], cutoff);
#endif

	  // left interval can be left out
	  feasRight = false;
	  oub [i] = ub [i] = fub;
	  chg_bds [i].setUpper (t_chg_bounds::CHANGED);
	  if (!(btCore (this, psi, cs, chg_bds, babInfo, true)))
	    goto end_procedure;
	}
      }

      // if both infeasible, fathom node
      if (!feasLeft && !feasRight && (fub < flb + COUENNE_EPS)) {
#ifdef DEBUG
	printf ("### prune all on infeasibility\n----------------------\n");
	for (int i=0; i<ncols; i++) printf ("   %2d %+20g %+20g  | %+20g\n",i,lb [i], ub [i], X [i]);
#endif
	goto end_procedure;
      }
    }

  retval = true;

end_procedure:

  // re-initialize to forget bound changes from previous tests
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
