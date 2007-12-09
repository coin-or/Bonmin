/*
 * Name:    checkNLP.cpp
 * Author:  Pietro Belotti
 * Purpose: temporary fix for checking NLP feasibility of incumbent integer solution
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"
#include "CouenneProblem.hpp"

// check if solution is MINLP feasible
bool CouenneProblem::checkNLP (const double *solution, const double obj) {

  // update variable array in evaluation structure
  //  expression::update (const_cast <double *> (solution), NULL, NULL);

  //  int nvars = extended ? nVars () : nOrig_;

  /*printf ("original NLP solution:\n");
  for (int i=0; i<nOrig_; i++)
    printf ("%4d: %10g\n", i, solution [i]);
    printf ("\n");*/

  CouNumber *sol = new CouNumber [nVars ()];
  CoinCopyN (solution, nOrig_, sol);

  /*printf ("copied NLP solution:\n");
  for (int i=0; i<nOrig_; i++)
    printf ("%4d: %10g\n", i, sol [i]);
    printf ("\n");*/

  getAuxs (sol);
  expression::update (sol, NULL, NULL);

  /*printf ("LP solution:\n");
  for (int i=0; i<nVars(); i++)
    printf ("%4d: %10g\n", i, sol [i]);
    printf ("\n");*/


  /*  CouNumber realobj = (*(p -> Obj (0) -> Body ())) ();
  if (fabs (realobj - obj) > COUENNE_EPS) {
    printf ("checkNLP: false objective function. %.3f != %.3f\n", realobj, obj);
    }*/

  // check (original and auxiliary) variables' integrality
  for (int i=0; i < nOrig_; i++) 

    if (variables_ [i] -> isInteger ()) {
      CouNumber val = expression::Variable (i);
      if (fabs (val - COUENNE_round (val)) > COUENNE_EPS) {
	Jnlst()->Printf(Ipopt::J_MOREDETAILED, J_PROBLEM,
			"checkNLP: integrality %d violated: %g [%g,%g]\n", 
			i, val, expression::Lbound (i), expression::Ubound (i));
	delete [] sol;
	return false;
      }
    }

  // check constraints
  if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM)) {
    for (int i=0; i < nCons (); i++) {

      CouenneConstraint *c = Con (i);

      CouNumber
	body = (*(c -> Body ())) (),
	lhs  = (*(c -> Lb   ())) (),
	rhs  = (*(c -> Ub   ())) ();

      if ((body > rhs + COUENNE_EPS) || 
	  (body < lhs - COUENNE_EPS)) {
	if (Jnlst()->ProduceOutput(Ipopt::J_WARNING, J_PROBLEM)) {
	  Jnlst()->Printf(Ipopt::J_WARNING, J_PROBLEM,
			  "Warning in checkNLP: constraint %d violated (lhs = %e body = %e rhs = %e): ", i, lhs, body, rhs);
	  c -> print ();
	}
	// We dont return anymore (return false;)
      }
    }
  }

  // check auxiliary variables
  //if (extended)
  /*
  if (0)
    for (int n = nVars (), i=0; i<n; i++) {

      int order = evalOrder (i);

      if (Var (order) -> Type () == AUX) {

	exprAux   *w   = dynamic_cast <exprAux *> (Var (order));
	CouNumber  
	  aux = expression::Variable (order),
	  img = (*(w -> Image ())) ();

	if (fabs (aux - img) > COUENNE_EPS) {
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
    }
  */

  delete [] sol;

  return true;
}
