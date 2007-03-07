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

inline CouNumber mymin (CouNumber a, CouNumber b) 
{return (a<b) ? a : b;} 

inline CouNumber mymax (CouNumber a, CouNumber b) 
{return (a>b) ? a : b;} 

/// Bound tightening for auxiliary variables

int CouenneCutGenerator::tightenBounds (const OsiSolverInterface &si, 
					char *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // Retrieve, from si, value and bounds of all variables, if not
  // firstcall, otherwise only those of the original ones Update
  // expression structure with x, l, u

  OsiSolverInterface *psi = const_cast <OsiSolverInterface *> (&si);

  CouNumber 
    *xc = const_cast <CouNumber *> (psi -> getColSolution ()),
    *lc = const_cast <CouNumber *> (psi -> getColLower    ()),
    *uc = const_cast <CouNumber *> (psi -> getColUpper    ());

  // update now all variables and bounds

  problem_ -> update (xc, lc, uc);

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. Here,
  // newly enforced branching rules may change dependent auxiliary
  // variables' bounds, in a recursive way. Hence we need to repeat
  // the propagation step as long as at least one bound is modified.

  bool found_one;

  do {

    found_one = false;

    int naux = problem_ -> nAuxs ();

    // check all auxiliary variables for changes in their upper,
    // lower bound, depending on the bound changes of the variables
    // they depend on

    for (register int i = problem_ -> nVars (), j=0; 
	 j < naux; j++) {

      CouNumber ll = (*(problem_ -> Aux (j) -> Lb ())) ();
      CouNumber uu = (*(problem_ -> Aux (j) -> Ub ())) ();

      //      printf ("x%3d: [%12.4f,%12.4f] -> [%12.4f,%12.4f] ", 
      //	      i+j, lc [i+j], uc [i+j], ll, uu);

      bool chg = false;

      // check if lower bound got higher    
      if (ll >= lc [i+j] + COUENNE_EPS) {
	//	psi -> setColLower (i+j, ll);
	//	printf ("l");
	lc [i+j] = ll;
	chg = true;
      }

      // check if upper bound got lower
      if (uu <= uc [i+j] - COUENNE_EPS) {
	//	psi -> setColUpper (i+j, uu);
	//	printf ("u");
	uc [i+j] = uu;
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

      expression::update (xc, lc, uc);
    }

  } while (found_one); // repeat as long as at least one bound changed

  // update again 
  problem_ -> update (xc, lc, uc);

  return nchg;
}
