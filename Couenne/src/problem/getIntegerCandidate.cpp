/*
 * Name:    getIntegerCandidate.cpp
 * Author:  Pietro Belotti
 * Purpose: generate integer NLP point Y starting from fractional 
 *          solution using bound tightening
 *
 * (C) Carnegie-Mellon University, 2007-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"

// lose patience after this many iterations of non-improving valid
// tightening (step 1)
#define VALID_ONLY_THRESHOLD 5 

/// generate integer NLP point xInt starting from fractional solution
/// xFrac, using (feasibility-based, i.e. cheap) bound tightening
///
/// return -1 if the problem is infeasible

int CouenneProblem::getIntegerCandidate (const double *xFrac, double *xInt, 
					 double *lb, double *ub) const {

  // return -1;

  jnlst_ -> Printf (J_MOREVECTOR, J_PROBLEM, "GetIntegerCandidate:---\n");

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
    objind   = Obj (0) -> Body () -> Index (),
    ncols    = nVars (), 
    retval   = 0,
    iter     = 0;

  double
    *olb   = new double [ncols],   *oub   = new double [ncols],  // outer bounds
    *llb   = new double [ncols],   *lub   = new double [ncols],  // new bounds when rounding down
    *rlb   = new double [ncols],   *rub   = new double [ncols],  // new                      up
    *dualL = new double [nOrig_],  *dualR = new double [nOrig_]; // lb[objind] per fix index/direction

  // copy in current bounds
  CoinCopyN (Lb (), ncols, olb);
  CoinCopyN (Ub (), ncols, oub);

  //fillFreeIntegers ();

  // start restricting around current integer box
  //for (int i=0; i<nVars(); i++) 
  for (int i=0; i<nOrig_; i++) 
    if ((Var (i) -> isInteger ()) &&    // integer, may fix if not dependent on other integers
	(true ||                        // TODO: remove and implement wave freeIntegers
	 (Var (i) -> Type () == VAR) || // if not aux
	 (freeIntegers_ [i]))) {        // 

      lb [i] = CoinMax (lb [i], floor (xFrac [i])); 
      ub [i] = CoinMin (ub [i], ceil  (xFrac [i]));
    }

  // for now save fractional point into integer point
  CoinCopyN (xFrac, nOrig_, xInt);

  domain_.push (nVars (), xInt, lb, ub);

  CoinCopyN (lb, nVars (), Lb ());
  CoinCopyN (ub, nVars (), Ub ());

  // translate current NLP point+bounds into higher-dimensional space
  initAuxs ();

  // TODO: re-copy first nOrig_ variables into xInt?
  //CoinCopyN (xFrac, nOrig_, xInt);

  // create fictitious change structure -- all initialized at UNCHANGED by constructor
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

  // structure to record fixed, non-fixed, and continuous variables
  typedef enum {UNFIXED, FIXED, CONTINUOUS} fixType;

  fixType *fixed = new fixType [nOrig_]; // integer variables that were fixed

  for (int i=0; i<nOrig_; i++) 
    fixed [i] = (Var (i) -> isInteger ()) ? UNFIXED : CONTINUOUS;

  bool changed;

  int validNonImproving = 0;

  // check if initAuxs() closed any bound; if so, fix variable
  for (int i=0; i<nOrig_; i++) 
    if (Var (i) -> isInteger () && (Lb (i) == Ub (i))) {
      X (i) = Lb (i);
      fixed [i] = FIXED;
    }

  do {

    if (jnlst_ -> ProduceOutput (Ipopt::J_MOREVECTOR, J_PROBLEM)) {

      printf ("===================================================\n");
      for (int i=0; i<nOrig_; i++)
	printf ("#### %4d: %d %c%c frac %20g  [%20g,%20g]\n", 
		i, fixed [i], 
		variables_ [i] -> isInteger () ? 'I' : ' ',
		(freeIntegers_ && freeIntegers_ [i]) ? 'F' : ' ',
		xFrac [i], Lb (i), Ub (i));
      printf ("---\n");
      for (int i=nOrig_; i<nVars (); i++)
	printf ("#### %4d:   %c%c frac %20g   [%20g,%20g]\n", 
		i, variables_ [i] -> isInteger () ? 'I' : ' ',
		(freeIntegers_ && freeIntegers_ [i]) ? 'F' : ' ',
		X (i), Lb (i), Ub (i));
      printf ("=================================================== %d %d\n", iter, nOrig_);
    }

    changed = false;

    for (int i = 0; i < nOrig_; i++)

      if (fixed [i] == UNFIXED) { // only check integer variables that have not been fixed yet

	// try rounding down ///////////////////////////////////////////////////////////////////////

	Lb (i) = Ub (i) = floor (xFrac [i]); 

	for (int j = 0; j<ncols; j++) {
	  f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
	  f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
	}

	f_chg [i].setLower (t_chg_bounds::CHANGED); 
	f_chg [i].setUpper (t_chg_bounds::CHANGED);

	bool feasLeft = btCore (f_chg); // true if feasible with fake bound

	dualL [i] = Lb (objind);

	// save new bounds 
	CoinCopyN (Lb (), ncols, llb);
	CoinCopyN (Ub (), ncols, lub);

	// restore initial situation
	CoinCopyN (olb, ncols, Lb ());
	CoinCopyN (oub, ncols, Ub ());

	// try rounding up ///////////////////////////////////////////////////////////////////////

	Lb (i) = Ub (i) = ceil (xFrac [i]); 

	for (int j = 0; j<ncols; j++) {
	  f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
	  f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
	}

	f_chg [i].setLower (t_chg_bounds::CHANGED); 
	f_chg [i].setUpper (t_chg_bounds::CHANGED);

	bool feasRight = btCore (f_chg); // true if feasible with fake bound

	dualR [i] = Lb (objind);

	// save new bounds
	CoinCopyN (Lb (), ncols, rlb);
	CoinCopyN (Ub (), ncols, rub);

	// restore initial situation
	CoinCopyN (olb, ncols, Lb ());
	CoinCopyN (oub, ncols, Ub ());

	//////////////////////////////////////////////////////////////////////////////////////////

	// Three cases:
	//
	// 1) if at least one is infeasible, set x_i to other
	//
	// 2) if both are infeasible, apply normal aggressive BT:
	//    2a) if both infeasible, node is infeasible
	//    2b) if both feasible, store index in free++ variables
	//    2c) if only one feasible, set rounded +/- 2
	//    ...
	//    2z) or probably simpler if return -1 to avoid calling Ipopt
	//
	// 3) if both feasible, choose one based on dual bound

	if (!feasLeft)

	  if (!feasRight) {
	    retval = -1; // case 2
	    break;
	  } else {

	    // ceil is feasible, floor is not.
	    fixed [i] = FIXED;
	    Lb (i) = Ub (i) = olb [i] = oub [i] = xInt [i] = ceil (xFrac [i]); 
	    //printf ("0 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
	    changed = true;
	    //printf ("+++ 1 %d\n", i);

	    // tighten bounds using r[lu]b
	    for (int j=0; j<ncols; j++) if (i != j) {

	      olb [j] = Lb (j) = CoinMax (Lb (j), rlb [j]);
	      oub [j] = Ub (j) = CoinMin (Ub (j), rub [j]);

	      if (Lb (j) > Ub (j) + COUENNE_EPS) {
		retval = -1;
		break;
	      }
	    }
	  }
	else if (!feasRight) {

	  // floor is feasible, ceil is not.

	  fixed [i] = FIXED;
	  Lb (i) = Ub (i) = olb [i] = oub [i] = xInt [i] = floor (xFrac [i]); 
	  //printf ("1 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
	  changed = true;
	  //printf ("+++ 2 %d\n", i);

	  // tighten bounds using l[lu]b
	  for (int j=0; j<ncols; j++) if (i != j) {

	    olb [j] = Lb (j) = CoinMax (Lb (j), llb [j]);
	    oub [j] = Ub (j) = CoinMin (Ub (j), lub [j]);

	    if (Lb (j) > Ub (j) + COUENNE_EPS) {
	      retval = -1;
	      break;
	    }
	  }
	} else { // case 3: tighten to smallest interval containing both [llb,lub] and [rlb,rub]

	  // tighten bounds using l[lu]b
	  for (int j=0; j<ncols; j++) {

	    olb [j] = Lb (j) = CoinMax (Lb (j), CoinMin (llb [j], rlb [j]));
	    oub [j] = Ub (j) = CoinMin (Ub (j), CoinMax (lub [j], rub [j]));

	    if (Lb (j) > Ub (j) + COUENNE_EPS) {
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
	    Lb (i) = Ub (i) = olb [i] = oub [i] = xInt[i] = floor (xFrac [i]);
	    fix = true;
	  }
	  else if (dualL [i] > dualR [i] + COUENNE_EPS) {
	    Lb (i) = Ub (i) = olb [i] = oub [i] = xInt[i] = ceil  (xFrac [i]);
	    fix = true;
	  }

	  if (fix) {
	    fixed [i] = FIXED;
	    //printf ("2 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
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

  if (retval >= 0)

    for (int i=0; i<nOrig_; i++) 
      if (fixed [i] == UNFIXED) {
	Lb (i) = Ub (i) = xInt [i] = 
	  ((dualL [i] < dualR [i] - COUENNE_EPS) ? floor (xFrac [i]) :
	   (dualL [i] > dualR [i] + COUENNE_EPS) ? ceil  (xFrac [i]) :
	   ((CoinDrand48 () < 0.5) ? floor (xFrac [i]) : ceil (xFrac [i])));

	//printf ("3 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
      }

  // save tightened bounds in NLP space
  CoinCopyN (Lb (), nOrig_, lb);
  CoinCopyN (Ub (), nOrig_, ub);

  // restore old bounds 
  //CoinCopyN (olb, ncols, Lb ());
  //CoinCopyN (oub, ncols, Ub ());

  for (int i=0; i<nOrig_; i++)
    if (fixed [i] == FIXED)
      lb [i] = ub [i] = xInt [i];

  if (jnlst_->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
    if (retval >= 0) {
      printf ("- retval %d ----------------------------------------------------------------\n", 
	      retval);
      for (int i=0; i<nOrig_; i++)
	printf ("#### %4d: %d %c frac %20g int %20g [%20g,%20g]\n", 
		i, fixed [i], variables_ [i] -> isInteger () ? 'I' : ' ',
		xFrac [i], xInt [i], lb [i], ub [i]);
    } else printf ("no good point was found\n");
  }

  delete [] f_chg;
  delete [] fixed;

  delete [] olb;   delete [] oub;
  delete [] dualL; delete [] dualR;
  delete [] llb;   delete [] lub;
  delete [] rlb;   delete [] rub;

  //printf (" done getintCand ------------------------------\n");

  domain_.pop ();

  jnlst_ -> Printf (J_MOREVECTOR, J_PROBLEM, "Done with GetIntegerCandidate\n");

  return retval;
}
