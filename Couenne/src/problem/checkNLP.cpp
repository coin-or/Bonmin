/*
 * Name:    checkNLP.cpp
 * Author:  Pietro Belotti
 * Purpose: check NLP feasibility of incumbent integer solution
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"

// check if solution is MINLP feasible
bool CouenneProblem::checkNLP (const double *solution, const double obj) {

  const int infeasible = 1;

  /*printf ("checking solution: [%g] ", obj);
  for (int i=0; i<nOrig_; i++)
    printf ("%.5f ", solution [i]);
    printf ("\n");*/

  CouNumber *sol = new CouNumber [nVars ()];

  // copy solution, evaluate the corresponding aux, and then replace
  // the original variables again for checking
  CoinCopyN (solution, nOrig_, sol);
  getAuxs (sol);
  CoinCopyN (solution, nOrig_, sol);

  // install NL solution candidate in evaluation structure
  domain_.push (nVars (), sol, domain_.lb (), domain_.ub ());

  /*printf ("checknlp: %d vars -------------------\n", domain_.current () -> Dimension ());
  for (int i=0; i<domain_.current () -> Dimension (); i++)
  printf ("%20g [%20g %20g]\n", domain_.x (i), domain_.lb (i), domain_.ub (i));*/

  CouNumber realobj = (*(Obj (0) -> Body ())) ();

  if (Obj (0) -> Sense () == MAXIMIZE)
    realobj = -realobj;

  bool retval = true;

  try {

    // check if objective corresponds

    if (fabs (realobj - obj) / (1 + fabs (realobj)) > feas_tolerance_) {

      Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
		      "checkNLP: false objective. %g != %g (diff. %g)\n", 
		      realobj, obj, realobj - obj);

      throw infeasible;
    }

    for (int i=0; i < nOrig_; i++) {

      CouNumber val = domain_.x (i);

      // check bounds

      if ((val > domain_.ub (i) + feas_tolerance_) ||
	  (val < domain_.lb (i) - feas_tolerance_)) {

	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: variable %d out of bounds: %.6f [%g,%g] (diff %g)\n", 
			i, val, domain_.lb (i), domain_.ub (i),
			CoinMax (fabs (val - domain_.lb (i)), 
				 fabs (val - domain_.ub (i))));
	throw infeasible;
      }

      // check (original and auxiliary) variables' integrality

      if (variables_ [i] -> isInteger () &&
	  (fabs (val - COUENNE_round (val)) > feas_tolerance_)) {

	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: integrality %d violated: %.6f [%g,%g]\n", 
			i, val, domain_.lb (i), domain_.ub (i));

	throw infeasible;
      }

      /*if (variables_ [i] -> Type () == AUX) {
	printf ("checking aux ");
	variables_ [i] -> print (); printf (" := ");
	variables_ [i] -> Image () -> print (); 
	printf (" --- %g = %g [%g]; {", 
		(*(variables_ [i])) (), 
		(*(variables_ [i] -> Image ())) (),
		(*(variables_ [i])) () -
		(*(variables_ [i] -> Image ())) ());
	for (int j=0; j<nVars (); j++)
	  printf ("%.6f ", (*(variables_ [j])) ());
	printf ("}\n");
	}*/

      CouNumber delta;

      // check if auxiliary has zero infeasibility

      if ((variables_ [i] -> Type () == AUX) && 
	  (fabs (delta = (*(variables_ [i])) () - 
		 (*(variables_ [i] -> Image ())) ()) > feas_tolerance_)) {

	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: auxiliarized %d violated (%g)\n", i, delta);

	throw infeasible;
      }
    }

    // check constraints

    if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM))

      for (int i=0; i < nCons (); i++) {

	CouenneConstraint *c = Con (i);

	CouNumber
	  body = (*(c -> Body ())) (),
	  lhs  = (*(c -> Lb   ())) (),
	  rhs  = (*(c -> Ub   ())) ();

	if ((body > rhs + feas_tolerance_ * (1 + CoinMax (fabs (body), fabs (rhs)))) || 
	    (body < lhs - feas_tolerance_ * (1 + CoinMax (fabs (body), fabs (rhs))))) {

	  if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM)) {

	    Jnlst()->Printf
	      (Ipopt::J_WARNING, J_PROBLEM,
	       "checkNLP: constraint %d violated (lhs=%+e body=%+e rhs=%+e, violation %g): ",
	       i, lhs, body, rhs, CoinMax (lhs-body, body-rhs));

	    c -> print ();
	  }

	  throw infeasible;
	  // We dont return anymore (return false;)
	}
      }
  }

  catch (int exception) {
    if (exception == infeasible) 
      retval = false;
  }

  delete [] sol;
  domain_.pop ();

  return retval;
}
