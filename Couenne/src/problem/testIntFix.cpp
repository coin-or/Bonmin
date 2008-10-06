/*
 * Name:    testIntFix.cpp
 * Author:  Pietro Belotti
 * Purpose: select rounding for integer variable based on tightening
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"

// test
int CouenneProblem::testIntFix (int index, 
				CouNumber xFrac, 
				enum fixType *fixed,
				CouNumber *xInt,
				CouNumber *dualL, CouNumber *dualR,
				CouNumber *olb,   CouNumber *oub,
				bool patient) const {
  int 
    ncols  = nVars (), 
    retval = 0,
    objind = Obj (0) -> Body () -> Index ();

  // create fictitious change structure -- all initialized at UNCHANGED by constructor
  t_chg_bounds *f_chg = new t_chg_bounds [ncols];

  double
    *llb = new double [ncols], *lub = new double [ncols],  // new bounds when rounding down
    *rlb = new double [ncols], *rub = new double [ncols];  //                          up

  // try rounding down ///////////////////////////////////////////////////////////////////////

  Lb (index) = Ub (index) = floor (xFrac); 

  for (int j = 0; j<ncols; j++) {
    f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
    f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
  }

  f_chg [index].setLower (t_chg_bounds::CHANGED); 
  f_chg [index].setUpper (t_chg_bounds::CHANGED);

  bool feasLeft = btCore (f_chg); // true if feasible with fake bound

  dualL [index] = Lb (objind);

  // save new bounds 
  CoinCopyN (Lb (), ncols, llb);
  CoinCopyN (Ub (), ncols, lub);

  // restore initial situation
  CoinCopyN (olb, ncols, Lb ());
  CoinCopyN (oub, ncols, Ub ());

  // try rounding up ///////////////////////////////////////////////////////////////////////

  Lb (index) = Ub (index) = ceil (xFrac); 

  for (int j = 0; j<ncols; j++) {
    f_chg [j].setLower (t_chg_bounds::UNCHANGED); 
    f_chg [j].setUpper (t_chg_bounds::UNCHANGED);
  }

  f_chg [index].setLower (t_chg_bounds::CHANGED); 
  f_chg [index].setUpper (t_chg_bounds::CHANGED);

  bool feasRight = btCore (f_chg); // true if feasible with fake bound

  dualR [index] = Lb (objind);

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

      jnlst_ -> Printf (J_MOREDETAILED, J_PROBLEM, 
			"test on %d -> Infeasible. ", index);
      retval = -1; // case 2

    } else {

      // ceil is feasible, floor is not.
      jnlst_ -> Printf (J_MOREDETAILED, J_PROBLEM, 
			"test on %d -> Right feasible, fix to %g. ", index, ceil (xFrac));

      fixed [index] = FIXED;
      Lb (index) = Ub (index) = olb [index] = oub [index] = xInt [index] = ceil (xFrac); 
      //printf ("0 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
      retval++;
      //printf ("+++ 1 %d\n", i);

      // tighten bounds using r[lu]b
      for (int j=0; j<ncols; j++) if (index != j) {

	olb [j] = Lb (j) = CoinMax (Lb (j), rlb [j]);
	oub [j] = Ub (j) = CoinMin (Ub (j), rub [j]);

	if (Lb (j) > Ub (j) + COUENNE_EPS)
	  retval = -1;
      }
    }
  else if (!feasRight) {

    // floor is feasible, ceil is not.
    jnlst_ -> Printf (J_MOREDETAILED, J_PROBLEM, 
		      "test on %d -> Left feasible, fix to %g. ", index, floor (xFrac));

    fixed [index] = FIXED;
    Lb (index) = Ub (index) = olb [index] = oub [index] = xInt [index] = floor (xFrac); 
    //printf ("1 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
    retval++;
    //printf ("+++ 2 %d\n", i);

    // tighten bounds using l[lu]b
    for (int j=0; j<ncols; j++) if (index != j) {

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

    if ((retval >= 0) && !patient) { // too much time spent here, just fix it based on dual bound

      fixed [index] = FIXED;

      Lb (index) = Ub (index) = olb [index] = oub [index] = xInt [index] = 
	((dualL [index] < dualR [index] - COUENNE_EPS) ? floor (xFrac) :
	 (dualL [index] > dualR [index] + COUENNE_EPS) ? ceil  (xFrac) :
	 ((CoinDrand48 () < 0.5) ? floor (xFrac) : ceil (xFrac)));
      
      jnlst_ -> Printf (J_MOREVECTOR, J_PROBLEM, 
			"test on %d -> Both feasible, lost patience, fixed to %g. ", 
			index, xInt [index]);

      //printf ("1 fixed %d [%g,%g,%g]\n", i, Lb (i), Ub (i), xInt [i]);
      retval++;
      //printf ("+++ 2 %d\n", i);
    } else if (retval >= 0) jnlst_ -> Printf (J_MOREVECTOR, J_PROBLEM, 
					      "test on %d -> Both feasible, skip this turn. ", index);
  }

  delete [] f_chg;

  delete [] llb; delete [] lub;
  delete [] rlb; delete [] rub;

  return retval;
}
