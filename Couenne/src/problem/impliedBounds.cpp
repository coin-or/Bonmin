/*
 * Name:    impliedBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: backward implied bound search
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CouenneProblem.h>

/// Bound tightening for auxiliary variables

int CouenneProblem::impliedBounds (char *chg_bds) const {

  int nchg = 0, //< number of bounds changed for propagation
      nvar = nVars ();

  /*for (int i=0; i < nVars () + nAuxs (); i++)
    printf ("x%3d: [%12.4f,%12.4f]\n", i, lb_ [i], ub_ [i]);*/

  for (int i=nAuxs (); i--;) {

    if (lb_ [nvar+i] > ub_ [nvar+i] + COUENNE_EPS) 
      return -1;

    //    if ((auxiliaries_ [i] -> Image () -> code () == COU_EXPRSUM) ||
    //	(auxiliaries_ [i] -> Image () -> code () == COU_EXPRGROUP))
    if (auxiliaries_ [i] -> Image () -> impliedBound (nvar+i, lb_, ub_, chg_bds) > COUENNE_EPS)
      nchg++;
  }

  /*  for (int i=0; i<nvar; i++) 
    printf (" [%g, %g]\n", 
	    expression::Lbound (i),
	    expression::Ubound (i));

  for (int i=0; i < nAuxs (); i++) {

    printf (" [%g, %g]", 
	    expression::Lbound (i+nvar),
	    expression::Ubound (i+nvar));

    auxiliaries_ [i] -> print (std::cout);
    printf (" := ");
    auxiliaries_ [i] -> Image () -> print (std::cout); fflush (stdout);
    printf ("\n");
    }*/

  return nchg;
}
