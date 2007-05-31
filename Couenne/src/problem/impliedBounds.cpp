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

int CouenneProblem::impliedBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0, //< number of bounds changed for propagation
      nvar = nVars ();

  /*printf ("=====================implied\n");
  for (int i=0; i < nVars () + nAuxs (); i++)
    if ((i < nVars ()) || (Aux (i-nVars ()) -> Multiplicity () > 0))
    printf ("x%d: [%g,%g]\n", i, lb_ [i], ub_ [i]);*/

  for (int i=nAuxs (); i--;) 

    //    if (i+nVars () != 79) // !!!
{

    //    for (int j=0; j<nAuxs () + nVars (); j++ )
    //      printf ("--- %d: [%.4f %.4f]\n", j, lb_ [j], ub_ [j]);

    if (lb_ [nvar+i] > ub_ [nvar+i] + COUENNE_EPS) {
      /*      printf ("w_%d has infeasible bounds [%g,%g]\n", 
	      i+nvar, lb_ [nvar+i], ub_ [nvar+i]);*/
      return -1;
    }

    //    if ((auxiliaries_ [i] -> Image () -> code () == COU_EXPRSUM) ||
    //	(auxiliaries_ [i] -> Image () -> code () == COU_EXPRGROUP))

    //    CouNumber l0 = lb_ [nvar+i], 
    //              u0 = ub_ [nvar+i];

   if (auxiliaries_ [i] -> Image () -> impliedBound (nvar+i, lb_, ub_, chg_bds) > COUENNE_EPS) {
     //      printf ("impli %2d [%g,%g] -> [%g,%g] ", nvar+i, l0, u0, lb_ [nvar+i], ub_ [nvar+i]);
     /*printf ("impli %2d ", nvar+i);

      auxiliaries_ [i] -> print (std::cout); printf (" := ");
      auxiliaries_ [i] -> Image () -> print (std::cout); printf ("\n");*/
      nchg++;
    }
  }
  /*
  for (int i=0; i<nvar; i++) 
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
    }
  */
  return nchg;
}
