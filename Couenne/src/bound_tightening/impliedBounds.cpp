/*
 * Name:    impliedBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: backward implied bound search
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"

//#define DEBUG

/// Bound tightening for auxiliary variables

int CouenneProblem::impliedBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

#ifdef DEBUG
  printf ("=====================implied\n");
  int j=0;
  for (int i=0; i < nVars (); i++) 
    if (variables_ [i] -> Multiplicity () > 0) {
      printf ("x_%03d [%+10g %+10g] ", i, 
	      expression::Lbound (i),
	      expression::Ubound (i));
      if (!(++j % 6)) printf ("\n");
    }
  if (j % 6) printf ("\n");
#endif

  for (int ii = nVars (); ii--;) {

    int i = numbering_ [ii];

    if (variables_ [i] -> Type () == AUX) {

      if (lb_ [i] > ub_ [i] + COUENNE_EPS) {
#ifdef DEBUG
	printf ("#### w_%d has infeasible bounds [%g,%g]\n", 
	i, lb_ [i], ub_ [i]);
#endif
	return -1;
      }

      //    if ((auxiliaries_ [i] -> Image () -> code () == COU_EXPRSUM) ||
      //	(auxiliaries_ [i] -> Image () -> code () == COU_EXPRGROUP))

#ifdef DEBUG
      CouNumber 
	l0 = lb_ [i], 
	u0 = ub_ [i];
#endif

      /*if (auxiliaries_ [i] -> Image () -> Argument () || 
	  auxiliaries_ [i] -> Image () -> ArgList  ()) {

	expression *arg = auxiliaries_ [i] -> Image () -> Argument ();
	if (!arg)   arg = auxiliaries_ [i] -> Image () -> ArgList  () [0];

	printf (":::: ");
	  arg -> print (std::cout);
	  if (arg -> Index () >= 0) {
	  int ind = arg -> Index ();
	  printf (" in [%g,%g]", 
	  expression::Lbound (ind), 
	  expression::Ubound (ind));
	  }
	  printf ("\n");
      }*/

      // TODO: also test if this expression, or any of its indep
      // variables, have changed. If not, skip

      if (variables_ [i] -> Image () -> impliedBound 
	  (variables_ [i] -> Index (), lb_, ub_, chg_bds)) {

#ifdef DEBUG
	printf ("impli %2d [%g,%g] -> [%g,%g]: ", i, l0, u0, lb_ [i], ub_ [i]);

	variables_ [i]             -> print (std::cout); printf (" := ");
	variables_ [i] -> Image () -> print (std::cout); printf ("\n");
#endif

	/*
	if (optimum_ && 
	    ((optimum_ [i+nvar] < lb_ [i+nvar] - COUENNE_EPS) ||
	     (optimum_ [i+nvar] > ub_ [i+nvar] + COUENNE_EPS)))
	  printf ("#### implied b_%d [%g,%g] cuts optimum %g: [%g --> %g, %g <-- %g]\n", 
		  i+nvar, expression::Lbound (i+nvar), expression::Ubound (i+nvar), 
		  optimum_ [i+nvar], l0, lb_ [i+nvar], ub_ [i+nvar], u0);

	//printf ("impli %2d ", nvar+i);

	if (auxiliaries_ [i] -> Image () -> Argument () || 
	    auxiliaries_ [i] -> Image () -> ArgList ()) {

	  expression *arg = auxiliaries_ [i] -> Image () -> Argument ();

	  if (!arg) {
	    for (int k=0; k < auxiliaries_ [i] -> Image () -> nArgs (); k++) {
	      arg =  auxiliaries_ [i] -> Image () -> ArgList () [k];
	      printf (" ");
	      arg -> print (std::cout);
	      if (arg -> Index () >= 0) {
		int ind = arg -> Index ();
		printf (" in [%g,%g]", 
			expression::Lbound (ind), 
			expression::Ubound (ind));
	      }	    
	    }
	  } else {
	    printf (" ");
	    arg -> print (std::cout);
	    if (arg -> Index () >= 0) {
	      int ind = arg -> Index ();
	      printf (" in [%g,%g]", 
		      expression::Lbound (ind), 
		      expression::Ubound (ind));
	    }
	  }
	} else printf (" [no args]");
	printf ("\n");
	*/

	nchg++;
      }
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
