/*
 * Name:    tightenBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: bound tightening for current linear relaxation
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>


/// Bound tightening for auxiliary variables

int CouenneProblem::tightenBounds (const OsiSolverInterface &si, 
				   char *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. Here,
  // newly enforced branching rules may change dependent auxiliary
  // variables' bounds, in a recursive way. Hence we need to repeat
  // the propagation step as long as at least one bound is modified.

  bool found_one;

  do {

    found_one = false;

    int naux = nAuxs ();

    // check all auxiliary variables for changes in their upper,
    // lower bound, depending on the bound changes of the variables
    // they depend on

    for (register int i = nVars (), j=0; 
	 j < naux; j++) {

      CouNumber ll = (*(Aux (j) -> Lb ())) ();
      CouNumber uu = (*(Aux (j) -> Ub ())) ();

      //      printf ("x%3d: [%12.4f,%12.4f] -> [%12.4f,%12.4f] ", 
      //	      i+j, lc [i+j], uc [i+j], ll, uu);

      bool chg = false;

      // check if lower bound got higher    
      if (ll >= lb_ [i+j] + COUENNE_EPS) {
	//	psi -> setColLower (i+j, ll);
	//	printf ("l");
	lb_ [i+j] = ll;
	chg = true;
      }

      // check if upper bound got lower
      if (uu <= ub_ [i+j] - COUENNE_EPS) {
	//	psi -> setColUpper (i+j, uu);
	//	printf ("u");
	ub_ [i+j] = uu;
	chg = true;
      }

      if (chg && chg_bds && !(chg_bds [i+j])) {
	nchg++;
	found_one = true;
	chg_bds [i+j] = 1; 
	//	printf (" ");
	//	problem_ -> Aux (j) -> Image () -> print (std::cout);
      }

      //      printf ("\n");

      // useless if assume expression::[lu]b etc already point to
      // problem::[lu]b
      expression::update (NULL, lb_, ub_);
    }

  } while (found_one); // repeat as long as at least one bound changed

  // update again 
  //  update (NULL, lb_, ub_);

  return nchg;
}
