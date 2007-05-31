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


/// Bound propagation for auxiliary variables

int CouenneProblem::tightenBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. 

  int naux = nAuxs ();

  // check all auxiliary variables for changes in their upper,
  // lower bound, depending on the bound changes of the variables
  // they depend on

  /*printf ("tighten========================\n");

  for (int i=0; i < nVars () + nAuxs (); i++) 
    if ((i < nVars ()) || (Aux (i-nVars ()) -> Multiplicity () > 0))
      printf ("x_%d [%g, %g]\n", i, 
	      expression::Lbound (i),
	      expression::Ubound (i));*/

  for (register int i = nVars (), j=0; 
       j < naux; j++) 

    if (Aux (j) -> Multiplicity () > 0) {

      //    for (int j=0; j<nAuxs () + nVars (); j++ )
      //      printf ("+++ %d: [%.4f %.4f]\n", j, lb_ [j], ub_ [j]);

      /*    printf ("w_%d [%g, %g] ", i+j,
	    expression::Lbound (i+j),
	    expression::Ubound (i+j));*/

      CouNumber ll = (*(Aux (j) -> Lb ())) (),
	uu = (*(Aux (j) -> Ub ())) ();

      //    printf (" ---> [%g, %g]\n ", ll, uu);

      /*auxiliaries_ [j] -> print (std::cout);
	printf (" := ");
	auxiliaries_ [j] -> Image () -> print (std::cout); fflush (stdout);
	printf ("\n");*/

      if (ll > uu + COUENNE_EPS) {
	/*
	printf ("w_%d has infeasible bounds [%g,%g]: ", i+j, ll, uu);
	Aux (j) -> Lb () -> print (std::cout); printf (" --- ");
	Aux (j) -> Ub () -> print (std::cout); printf ("\n");
	*/
	return -1; // declare this node infeasible
      }

      //      bool chg = false;

      // check if lower bound got higher    
      if ((ll > - COUENNE_INFINITY) && (ll >= lb_ [i+j] + COUENNE_EPS)) {

	/*printf ("update lbound %d: %.10e >= %.10e + %.12e\n", 
	  i+j, ll, lb_ [i+j], COUENNE_EPS);*/

	/*printf ("update lbound %d: %g >= %g\n", 
	  i+j, ll, lb_ [i+j]);*/

	/*printf ("propa %2d [%g,(%g)] -> [%g,(%g)] (%g)", 
		i+j, lb_ [i+j], ub_ [i+j], ll, uu, lb_ [i+j] - ll);
	Aux (j)             -> print (std::cout); printf (" := ");
	Aux (j) -> Image () -> print (std::cout); printf ("\n");*/

	lb_ [i+j] = ll;
	chg_bds [i+j].lower = CHANGED;
	nchg++;
      }

      // check if upper bound got lower
      if ((uu < COUENNE_INFINITY) && (uu <= ub_ [i+j] - COUENNE_EPS)) {

	/*printf ("update ubound %d: %.10e <= %.10e - %.12e (%.12e)\n", 
	  i+j, uu, ub_ [i+j], COUENNE_EPS, uu - ub_ [i+j]);*/
	/*printf ("update ubound %d: %g >= %g\n", 
	  i+j, uu, ub_ [i+j]);*/

	/*printf ("propa %2d [(%g),%g] -> [(%g),%g] (%g)", 
		i+j, lb_ [i+j], ub_ [i+j], ll, uu, ub_ [i+j] - uu);
	Aux (j)             -> print (std::cout); printf (" := ");
	Aux (j) -> Image () -> print (std::cout); printf ("\n");*/

	ub_ [i+j] = uu;
	chg_bds [i+j].upper = CHANGED;
	nchg++;
      }

      // useless if assume expression::[lu]b_ etc already point to
      // problem::[lu]b_
      expression::update (NULL, lb_, ub_);
    }

  return nchg;
}
