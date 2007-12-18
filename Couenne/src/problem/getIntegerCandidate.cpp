/*
 * Name:    getIntegerCandidate.cpp
 * Author:  Pietro Belotti
 * Purpose: generate integer NLP point Y starting from fractional 
 *          solution using bound tightening
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"

// lose patience after this many iterations of non-improving valid
// tightening (step 1)
#define VALID_ONLY_THRESHOLD 5 

/// generate integer NLP point Y starting from fractional solution
/// using bound tightening
///
/// return -1 if the problem is infeasible

int CouenneProblem::getIntegerCandidate (const double *xFrac, double *xInt, 
					 double *lb, double *ub) {

  // A more sophisticated rounding procedure
  //
  // Call aggressive tightening for all integer variables, setting
  // each first at floor(x_i) and then at ceil(x_i):
  //
  // 1) if at least one is infeasible, set x_i to other
  //
  // 2) if both are infeasible, apply normal aggressive BT:
  //    2a) if both infeasible, node is infeasible
  //    2b) if both feasible, store index in free++ variables
  //    2c) if only one feasible, set rounded +/- 2
  //
  // 3) if both feasible, store index in free variables
  //
  // with free and free++ variables, look for integer point closest to
  // xcurrent.

  int
    objind   = Obj (0) -> Body  () -> Index (),
    ncols    = nVars (), 
    retval   = 0,
    iter     = 0;

  double
    *olb   = new double [ncols],   *oub   = new double [ncols],  // saved bounds
    *llb   = new double [ncols],   *lub   = new double [ncols],  // bounds when rounding down
    *rlb   = new double [ncols],   *rub   = new double [ncols],  //                      up
    *dualL = new double [nOrig_],  *dualR = new double [nOrig_]; // lb[objind] per fix index/direction

  // save old bounds
  CoinCopyN (lb_, ncols, olb);
  CoinCopyN (ub_, ncols, oub);

  // for now save fractional point into integer point
  CoinCopyN (xFrac, nOrig_, xInt);

  // translate current NLP point+bounds into higher-dimensional space
  initAuxs (xInt, lb, ub);

  // create fictitious change structure -- all initialized at UNCHANGED by constructor
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

  // structure to record fixed, non-fixed, and continuous variables
  typedef enum {UNFIXED, FIXED, CONTINUOUS} fixType;
  fixType *fixed = new fixType [nOrig_]; // integer variables that were fixed
  for (int i=0; i<nOrig_; i++) fixed [i] = (Var (i) -> isInteger ()) ? UNFIXED : CONTINUOUS;

  bool changed;

  /*printf ("========================================================================\n");
  printf ("========================================================================\n");
  for (int i=0; i<nOrig_; i++)
    printf ("#### %4d: %d %c frac %20g                          [%20g,%20g]\n", 
	    i, fixed [i], variables_ [i] -> isInteger () ? 'I' : ' ',
	    xFrac [i], lb_ [i], ub_ [i]);*/

  int validNonImproving = 0;

  do {

    changed = false;

    // TODO: use explicit checkNLP for solutions

    //printf ("=================================================== %d %d\n", iter, nOrig_);

    for (int i = 0; i < nOrig_; i++)

      if (fixed [i] == UNFIXED) { // only check integer variables that have not been fixed yet

	// try rounding down ///////////////////////////////////////////////////////////////////////

	lb_ [i] = ub_ [i] = floor (xFrac [i]); 

	for (int j = 0; j<ncols; j++) {
	  f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
	  f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
	}

	f_chg [i].setLower (t_chg_bounds::CHANGED); 
	f_chg [i].setUpper (t_chg_bounds::CHANGED);

	bool feasLeft = btCore (f_chg); // true if feasible with fake bound

	dualL [i] = lb_ [objind];

	// save new bounds 
	CoinCopyN (lb_, ncols, llb);
	CoinCopyN (ub_, ncols, lub);

	// restore initial situation
	CoinCopyN (olb, ncols, lb_);
	CoinCopyN (oub, ncols, ub_);

	// try rounding up ///////////////////////////////////////////////////////////////////////

	lb_ [i] = ub_ [i] = ceil (xFrac [i]); 

	for (int j = 0; j<ncols; j++) {
	  f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
	  f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
	}

	f_chg [i].setLower (t_chg_bounds::CHANGED); 
	f_chg [i].setUpper (t_chg_bounds::CHANGED);

	bool feasRight = btCore (f_chg); // true if feasible with fake bound

	dualR [i] = lb_ [objind];

	// save new bounds
	CoinCopyN (lb_, ncols, rlb);
	CoinCopyN (ub_, ncols, rub);

	// restore initial situation
	CoinCopyN (olb, ncols, lb_);
	CoinCopyN (oub, ncols, ub_);

	// Three cases:
	//
	// 1) if at least one is infeasible, set x_i to other
	//
	// 2) if both are infeasible, apply normal aggressive BT:
	//    2a) if both infeasible, node is infeasible
	//    2b) if both feasible, store index in free++ variables
	//    2c) if only one feasible, set rounded +/- 2
	//    ...
	//    2z) or probably simpler if return -1 to tell our hero not to
	//    call Ipopt
	//
	// 3) if both feasible, choose one based on dual bound

	if (!feasLeft)

	  if (!feasRight) {
	    retval = -1; // case 2
	    break;
	  } else {

	    // ceil is feasible, floor is not.
	    fixed [i] = FIXED;
	    lb_ [i] = ub_ [i] = xInt [i] = ceil (xFrac [i]); 
	    changed = true;
	    //printf ("+++ 1 %d\n", i);

	    // tighten bounds using r[lu]b
	    for (int j=0; j<ncols; j++) if (i != j) {

	      lb_ [j] = CoinMax (lb_ [j], rlb [j]);
	      ub_ [j] = CoinMin (ub_ [j], rub [j]);

	      if (lb_ [j] > ub_ [j] + COUENNE_EPS) {
		retval = -1;
		break;
	      }
	    }
	  }
	else if (!feasRight) {

	  // floor is feasible, ceil is not.

	  fixed [i] = FIXED;
	  lb_ [i] = ub_ [i] = xInt [i] = floor (xFrac [i]); 
	  changed = true;
	  //printf ("+++ 2 %d\n", i);

	  // tighten bounds using l[lu]b
	  for (int j=0; j<ncols; j++) if (i != j) {

	    lb_ [j] = CoinMax (lb_ [j], llb [j]);
	    ub_ [j] = CoinMin (ub_ [j], lub [j]);

	    if (lb_ [j] > ub_ [j] + COUENNE_EPS) {
	      retval = -1;
	      break;
	    }
	  }
	} else { // case 3: tighten to smallest interval containing both [llb,lub] and [rlb,rub]

	  // tighten bounds using l[lu]b
	  for (int j=0; j<ncols; j++) {

	    lb_ [j] = CoinMax (lb_ [j], CoinMin (llb [j], rlb [j]));
	    ub_ [j] = CoinMin (ub_ [j], CoinMax (lub [j], rub [j]));

	    if (lb_ [j] > ub_ [j] + COUENNE_EPS) {
	      retval = -1;
	      break;
	    }
	  }
	}
      }

    if (retval == -1) break;

    if (!changed) {

      validNonImproving++;

      for (int i=0; i<nOrig_; i++) 
	if (fixed [i] == UNFIXED) {

	  bool fix = false;

	  if      (dualL [i] < dualR [i] - COUENNE_EPS) {
	    lb_[i] = ub_[i] = xInt[i] = floor (xFrac [i]);
	    fix = true;
	  }
	  else if (dualL [i] > dualR [i] + COUENNE_EPS) {
	    lb_[i] = ub_[i] = xInt[i] = ceil  (xFrac [i]);
	    fix = true;
	  }

	  if (fix) {
	    fixed [i] = FIXED;
	    //printf ("+++ 3 --- %d %d %d\n", i, validNonImproving, VALID_ONLY_THRESHOLD);
	    changed = true;
	    if (validNonImproving < VALID_ONLY_THRESHOLD) break;
	    //else printf ("fix more... ");
	  }
	}
    }

  } while ((iter++ < nOrig_) && changed && (validNonImproving < 2 * VALID_ONLY_THRESHOLD));

  // phase 2: the fixable has been fixed, now the unfixed integer
  // variables have to be fixed too, e.g. according to dual bound

  for (int i=0; i<nOrig_; i++) 
    if (fixed [i] == UNFIXED)
      lb_ [i] = ub_ [i] = xInt [i] = 
	((dualL [i] < dualR [i] - COUENNE_EPS) ? floor :
	 (dualL [i] > dualR [i] + COUENNE_EPS) ? ceil  :
	 ((CoinDrand48 () < 0.5) ? floor : ceil)) (xFrac [i]);

  // save tightened bounds in NLP space
  CoinCopyN (lb_, nOrig_, lb);
  CoinCopyN (ub_, nOrig_, ub);

  // restore old bounds 
  CoinCopyN (olb, ncols, lb_);
  CoinCopyN (oub, ncols, ub_);

  /*printf ("------------------------------------------------------------------------\n");
  for (int i=0; i<nOrig_; i++)
    printf ("#### %4d: %d %c frac %20g int %20g [%20g,%20g]\n", 
	    i, fixed [i], variables_ [i] -> isInteger () ? 'I' : ' ',
	    xFrac [i], xInt [i], lb_ [i], ub_ [i]);
	    printf ("========================================================================\n");*/

  delete [] olb;
  delete [] oub;
  delete [] f_chg;

  delete [] dualL;
  delete [] dualR;
  delete [] llb;
  delete [] lub;
  delete [] rlb;
  delete [] rub;

  delete [] fixed;

  //printf (" done getintCand ------------------------------\n\n\n");

  return retval;
}
