/*
 * Name:    checkNLP.cpp
 * Author:  Pietro Belotti
 * Purpose: temporary fix for checking NLP feasibility of incumbent integer solution
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */


#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

// check if solution is MINLP feasible
bool checkNLP (CglCutGenerator *g, const double *solution, double &obj) {

  // first cut generator (if the only one) is a CouenneCutGenerator,
  // which has NLP symbolic information. Use that to check NLP feasibility
  CouenneCutGenerator *cg = dynamic_cast <CouenneCutGenerator *> (g);
  if (!cg) return false;

  CouenneProblem *p = cg -> Problem ();
  if (!p) return false;

  // update variable array in evaluation structure
  expression::update (const_cast <double *> (solution), NULL, NULL);

  /*  CouNumber realobj = (*(p -> Obj (0) -> Body ())) ();
  if (fabs (realobj - obj) > COUENNE_EPS) {
    printf ("checkNLP: false objective function. %.3f != %.3f\n", realobj, obj);
    }*/

  // check (original and auxiliary) variables' integrality
  for (int i=0; i < p -> nVars (); i++) 

    if (p -> Var (i) -> isInteger ()) {
      CouNumber val = expression::Variable (i);
      if (fabs (val - COUENNE_round (val)) > COUENNE_EPS)
	return false;
    }

  // check constraints

  for (int i=0; i < p -> nCons (); i++) {

    CouenneConstraint *c = p -> Con (i);

    CouNumber body = (*(c -> Body ())) (),
              lhs  = (*(c -> Lb   ())) (),
              rhs  = (*(c -> Ub   ())) ();

    if ((body > rhs + COUENNE_EPS) || 
	(body < lhs - COUENNE_EPS))
      return false;
  }

  // check auxiliary variables

  for (int n = p -> nVars (), i=0; i<n; i++) {

    int order = p -> evalOrder (i);

    if (p -> Var (order) -> Type () == AUX) {

      exprAux   *w   = dynamic_cast <exprAux *> (p -> Var (i));
      CouNumber  aux = expression::Variable (i),
	         img = (*(w -> Image ())) ();

      if (fabs (aux - img) > COUENNE_EPS) {

	//	printf ("checkNLP: Auxiliary different from its expression\n");
	return false;
      } 
    }
  }

  return true;
}
