/*
 * Name:    tightenBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: bound tightening for current linear relaxation
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

//#define DEBUG

/// Bound propagation for auxiliary variables

int CouenneProblem::tightenBounds (t_chg_bounds *chg_bds) const {

  int nchg = 0; //< number of bounds changed for propagation

  // update bounding box (which may depend on the original
  // variables' box) for the variables whose bound is looser. 

  // check all auxiliary variables for changes in their upper,
  // lower bound, depending on the bound changes of the variables
  // they depend on

  if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
    // ToDo: Pipe all output through journalist
    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
		    "  tighten========================\n  ");
    int j=0;
    for (int i=0; i < nVars (); i++) 
      if (variables_ [i] -> Multiplicity () > 0) {
	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"x_%03d [%+10g %+10g] ", i, 
			domain_. lb (i),
			domain_. ub (i));
	if (!(++j % 6)) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n  ");
    }
    if (j % 6) Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
  }

  // Must go -- but change it first!
  //expression::update (NULL, Lb (), Ub ());

  for (register int ii = 0, j = nVars (); j--; ii++) {

    int i = numbering_ [ii];

    // early test to avoid a loop

    if (Lb (i) > Ub (i) + COUENNE_EPS) {

      if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			"pre-check: w_%d has infeasible bounds [%g,%g]. ", i, Lb (i), Ub (i));
	Var (i) -> Lb () -> print (std::cout);
	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," --- ");
	Var (i) -> Ub () -> print (std::cout);
	Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
      }

      return -1; // declare this node infeasible
    }

    //printf ("w_%d [%g, %g] ", i,
    //domain_. lb (i),
    //domain_. ub (i));

    if ((Var (i) -> Multiplicity () > 0) &&
	(Var (i) -> Type         () == AUX) 
	// TODO: also test if any indep variable of this expression
	// have changed. If not, skip
	) {

      CouNumber ll, uu; 
      //ll = (*(variables_ [i] -> Lb ())) (),
      //uu = (*(variables_ [i] -> Ub ())) ();

      variables_ [i] -> Image () -> getBounds (ll, uu);

      /*printf (" ---> [%g, %g] ", ll, uu);
      variables_ [i] -> print (std::cout);
      printf (" := ");
      variables_ [i] -> Image () -> print (std::cout); fflush (stdout);
      printf (" [");
      variables_ [i] -> Lb () -> print (std::cout); fflush (stdout);
      printf (" , ");
      variables_ [i] -> Ub () -> print (std::cout); fflush (stdout);
      printf ("]\n");*/

      if (ll > uu + COUENNE_EPS) {

	if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			  "w_%d has infeasible bounds [%g,%g]: ", i, ll, uu);
	  Var (i) -> Lb () -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," --- ");
	  Var (i) -> Ub () -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
	}

	return -1; // declare this node infeasible
      }

      // check if lower bound got higher
      if ((ll > - COUENNE_INFINITY) && 
	  (ll >= Lb (i) + COUENNE_EPS) &&
	  ((fabs (ll)        < COUENNE_EPS) || 
	   (fabs (Lb (i)) < COUENNE_EPS) ||
	   (fabs (ll / (Lb (i)) - 1) > COUENNE_EPS)) ) {

	/*if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {

	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			  "propa %2d [%g,(%g)] -> [%g,(%g)] (%g) ", 
			  i, Lb (i), Ub (i), ll, uu, Lb (i) - ll);
	  Var (i)             -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," := ");
	  Var (i) -> Image () -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");

	  if (optimum_ && 
	      (optimum_ [i] >= Lb (i)) && 
	      (optimum_ [i] <= ll - COUENNE_EPS)) {

	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			    "#### propagating l_%d cuts optimum: [%g --> %g -X-> %g] :: ", 
			    i+j, Lb (i), optimum_ [i], ll);
	    Var (i) -> Lb () -> print (std::cout);
	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," --- ");
	    Var (i) -> Ub () -> print (std::cout);
	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
	  }
	  }*/

	Lb (i) = ll;
	chg_bds [i].setLower(t_chg_bounds::CHANGED);
	nchg++;
      }

      // check if upper bound got lower
      if ((uu < COUENNE_INFINITY) && 
	  (uu <= Ub (i) - COUENNE_EPS) &&
	  ((fabs (uu)      < COUENNE_EPS) || 
	   (fabs (Ub (i)) < COUENNE_EPS) ||
	   (fabs (uu / (Ub (i)) - 1) > COUENNE_EPS)) ) {
	//      if ((uu < COUENNE_INFINITY) && (uu <= ub_ [i+j] - COUENNE_EPS)) {

	/*printf ("update ubound %d: %.10e <= %.10e - %.12e (%.12e)\n", 
	  i+j, uu, ub_ [i+j], COUENNE_EPS, uu - ub_ [i+j]);*/
	/*printf ("update ubound %d: %g >= %g\n", 
	  i+j, uu, ub_ [i+j]);*/

	if (Jnlst()->ProduceOutput(J_VECTOR, J_BOUNDTIGHTENING)) {
	  /*Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			  "prop %2d [(%g),%g] -> [(%g),%g] (%g) ", 
			  i, Lb (i), Ub (i), ll, uu, Ub (i) - uu);
	  Var (i)             -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," := ");
	  Var (i) -> Image () -> print (std::cout);
	  Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");*/

	  if (optimum_ && 
	      (optimum_ [i] <= Ub (i)) && 
	      (optimum_ [i] >= uu + COUENNE_EPS)) {

	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,
			    "Couenne: propagating u_%d cuts optimum: [%g <-X- %g <-- %g] :: ", 
			    i, uu, optimum_ [i], Ub (i));
	    Var (i) -> Lb () -> print (std::cout);
	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING," --- ");
	    Var (i) -> Ub () -> print (std::cout);
	    Jnlst()->Printf(J_VECTOR, J_BOUNDTIGHTENING,"\n");
	  }
	}

	Ub (i) = uu;
	chg_bds [i].setUpper(t_chg_bounds::CHANGED);
	nchg++;
      }

      // useless if assume expression::[lu]b_ etc already point to
      // problem::[lu]b_
    }
  }

  if (nchg)
    Jnlst () -> Printf (J_VECTOR, J_BOUNDTIGHTENING,
			"  forward tightening %d changes\n", nchg);

  return nchg;
}
