/*
 * Name:    checkNLP.cpp
 * Author:  Pietro Belotti
 * Purpose: temporary fix for checking NLP feasibility of incumbent integer solution
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */


#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

bool checkNLP (CglCutGenerator *g, const double *solution, double &obj) {

  // PBelotti: check NLP feasibility

  // first cut generator (if the only one) is a CouenneCutGenerator,
  // which has NLP symbolic information. Use that to check NLP feasibility
    
  CouenneProblem *p = dynamic_cast <CouenneCutGenerator *> (g) -> Problem ();

  // update variable array in evaluation structure

  expression::update (const_cast <double *> (solution), NULL, NULL);

  // check constraints

  CouNumber realobj = (*(p -> Obj (0) -> Body ())) ();

  /*  printf ("objective: %.4f = %.4f = ", obj, realobj); 
  p -> Obj (0) -> Body () -> print (std::cout);
  printf (" (variable %d)\n", p -> Obj (0) -> Body () -> Index ());*/

  if (fabs (realobj - obj) > COUENNE_EPS) {
    obj = realobj;
  }

  for (int i=0; i<p->nNLCons (); i++) {

    CouenneConstraint *c = p -> NLCon (i);

    CouNumber body = (*(c -> Body ())) (),
              lhs  = (*(c -> Lb   ())) (),
              rhs  = (*(c -> Ub   ())) ();

    //    printf ("is %.2f <= %.2f <= %.2f? ", lhs, body, rhs);

    if ((body > rhs + COUENNE_EPS) || 
	(body < lhs - COUENNE_EPS)) {

      //      printf ("############### NO!\n");
      return false;
    }
    //    else printf ("yes\n");
  }

  // check auxiliary variables

  for (int i=0; i<p->nAuxs (); i++) {

    exprAux *w = p -> Aux (i);

    CouNumber aux = expression::Variable (i+p->nVars ()),
              img = (*(w -> Image ())) ();

    //    printf ("is w_%d = %.2f really equal to %.2f =", i+p->nVars(), aux, img);
    //    w -> Image () -> print (std::cout);
    //    printf ("? ");

    if (fabs (aux - img) > COUENNE_EPS) {

      //      printf ("######### NO!\n");
      return false;
    } 
    //    else printf ("yes (%.15f)\n", aux-img);
  }

  return true;
}
  
