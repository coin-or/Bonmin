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

#define FEAS_TOL 1e-5

// exception
#define INFEASIBLE 1

// check if solution is MINLP feasible
bool CouenneProblem::checkNLP (const double *solution, const double obj) {

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

  domain_.push (nVars (), sol, domain_.lb (), domain_.ub ());
  //domain_.current () -> resize (nVars ());

  /*printf ("checknlp: %d vars -------------------\n", domain_.current () -> Dimension ());
  for (int i=0; i<domain_.current () -> Dimension (); i++)
  printf ("%20g [%20g %20g]\n", domain_.x (i), domain_.lb (i), domain_.ub (i));*/

  CouNumber realobj = (*(Obj (0) -> Body ())) ();

  if (Obj (0) -> Sense () == MAXIMIZE)
    realobj = -realobj;

  try {

    if (fabs (realobj - obj) > FEAS_TOL) {
      Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
		      "checkNLP: warning, false objective. %.3f != %.3f (%g)\n", 
		      realobj, obj, realobj - obj);
    }

    // check (original and auxiliary) variables' integrality
    for (int i=0; i < nOrig_; i++) {

      CouNumber val = domain_.x (i);//expression::Variable (i);

      if ((val > domain_.ub (i) + FEAS_TOL) ||
	  (val < domain_.lb (i) - FEAS_TOL)) {

	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: variable %d out of bounds: %.6f [%g,%g] (diff %g)\n", 
			i, val, domain_.lb (i), domain_.ub (i),
			CoinMax (fabs (val - domain_.lb (i)), 
				 fabs (val - domain_.ub (i))));
	throw INFEASIBLE;
      }

      if (variables_ [i] -> isInteger () &&
	  (fabs (val - COUENNE_round (val)) > FEAS_TOL)) {
	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: integrality %d violated: %.6f [%g,%g]\n", 
			i, val, domain_.lb (i), domain_.ub (i));
	throw INFEASIBLE;
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

      if ((variables_ [i] -> Type () == AUX) && 
	  (fabs (delta = (*(variables_ [i])) () - (*(variables_ [i] -> Image ())) ()) > FEAS_TOL)) {
	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: auxiliarized %d violated (%g)\n", i, delta);
	throw INFEASIBLE;
      }
    }

    // check constraints
    if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM))

      for (int i=0; i < nCons (); i++) {

	CouenneConstraint *c = Con (i);

	/*printf ("checking constraint ");
	  c -> print ();*/

	CouNumber
	  body = (*(c -> Body ())) (),
	  lhs  = (*(c -> Lb   ())) (),
	  rhs  = (*(c -> Ub   ())) ();

	if ((body > rhs + FEAS_TOL) || 
	    (body < lhs - FEAS_TOL)) {

	  if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM)) {

	    Jnlst()->Printf
	      (Ipopt::J_WARNING, J_PROBLEM,
	       "Warning in checkNLP: constraint %d violated (lhs=%+e body=%+e rhs=%+e, violation %g): ",
	       i, lhs, body, rhs, CoinMax (lhs-body, body-rhs));

	    c -> print ();
	  }
	  // We dont return anymore (return false;)
	}
      }

    // check auxiliary variables
    //if (extended)

    /*for (int n = nVars (), i=0; i<n; i++) {

      int order = evalOrder (i);

      if (Var (order) -> Type () == AUX) {

	exprAux   *w   = dynamic_cast <exprAux *> (Var (order));
	CouNumber  
	  aux = expression::Variable (order),
	  img = (*(w -> Image ())) ();

	if (fabs (aux - img) > FEAS_TOL) {
	  if (Jnlst()->ProduceOutput(Ipopt::J_MOREDETAILED, J_PROBLEM)) {
	    Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			    "auxiliary %d infeasible: ", order);
	    w -> print ();
	    Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM," := ");
	    w -> Image () -> print ();
	    Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,"\n");
	  }
	  return false;
	}
      }
    }*/

    delete [] sol;
    domain_.pop ();
  }

  catch (int exception) {
    if (exception == INFEASIBLE) 
      return false;
  }

  return true;
}
