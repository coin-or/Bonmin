/*
 * Name:    standardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize all expressions in a problem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include <CoinHelperFunctions.hpp>
#include <CoinTime.hpp>

#include <CouenneTypes.h>

#include <expression.hpp>
#include <exprClone.hpp>

#include <CouenneProblem.hpp>
#include <CouenneProblemElem.hpp>
#include <depGraph.hpp>


/// standardize (nonlinear) common expressions, objectives, and constraints

void CouenneProblem::standardize () {

  graph_ = new DepGraph;

  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); i++)
    graph_ -> insert (*i);

  // allocate space in auxiliaries_ from commonexprs_

  /*for (std::vector <expression *>::iterator i = commonexprs_.begin (); 
       i != commonexprs_.end (); i++) 
       auxiliaries_. push_back (NULL);*/

  // standardize initial aux variables (aka defined variables, common
  // expression)
  /*for (int i=0; i < commonexprs_ . size (); i++) {

    exprAux *aux = addAuxiliary (commonexprs_ [i]);
    //exprAux *aux = auxiliaries_ [i];

    exprAux *w = new exprAux (symbolic, 
			    variables_ . size () + auxiliaries_ . size (), 
			    1 + symbolic -> rank (this));

    printf ("////////////// now attempting to standardize defVar "); fflush (stdout);
    aux -> print ();
    printf (" := "); fflush (stdout);
    aux -> Image () -> print (); 
    printf ("\n ----> "); fflush (stdout);

    exprAux *naux = aux -> Image () -> standardize (this);

    if (naux) {
      printf ("done: "); fflush (stdout);
      naux -> print (); printf ("\n");
      //printf (" := "); fflush (stdout);
      //naux -> Image () -> print (); printf ("\n..."); fflush (stdout);
    } else if (aux) {
      aux -> print ();
      printf (" := "); fflush (stdout);
      aux -> Image () -> print (); printf ("\n");
    } else printf ("[n]aux NULL!\n");
    //    if (aux) 
    //  (*i) -> Body (new exprClone (aux));
    }*/

  // standardize objectives
  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);
    if (aux) 
      (*i) -> Body (new exprClone (aux));
  }

  // commuted_ is an array with a flag for each original variable,
  // which is true at position i if initially original variable x_i
  // went auxiliary

  commuted_ = new bool [nVars ()];
  for (int i = nVars (); i--;)
    *commuted_++ = false;
  commuted_ -= nVars ();

  std::vector <CouenneConstraint *> con2;


  // standardize constraints
  for (int i=0, ncon = constraints_.size (); ncon--; i++) {

    /*printf ("prima  --------------------------------------\n");
    print ();
    printf ("--------------------------------------\n");*/

    exprAux *aux = constraints_ [i] -> standardize (this);

    /*printf ("dopo--------------------------------------\n");
      print ();*/

    //Back to normal
    /*if (!aux) {
      delete constraints_ [i];
      constraints_ [i] = NULL;
      continue;
      }*/

    if (aux) { // save if standardized
      constraints_ [i] -> Body (new exprClone (aux));
      con2.push_back (constraints_ [i]);
    }

    // now aux is an auxiliary variable. If it is linear, three cases:
    //
    // 1) it is an expr{Sum,Group,Sub}
    //
    //  1a) it is an equality constraint. Look for a proper element to
    //      "recreate" an auxiliary variable from within the linear
    //      expression, and discard the one just created or simply
    //      decrease its multiplicity
    //
    //  1b) it is an inequality, just replace the body of the
    //      constraint
    //
    // 2) it is a variable itself, which means the original constraint
    //    was of the form aw >= k or w = h. In both cases, the auxiliary
    //    variable can be removed and replaced with its image
    //
    // 3) it is nonlinear: the good old case, do nothing.

    /*printf ("=== "); fflush (stdout); 
    aux -> print (); printf (" := "); fflush (stdout);
    aux -> Image () -> print (); printf ("\n");*/
  }

  //  constraints_. erase (constraints_.begin (), constraints_.end ());

  constraints_ = con2;

  //print ();

  delete auxSet_;

  //printf ("ntotvars = %d\n", nVars());

  int nTotVar = nVars ();// + nAuxs ();

  // reallocate space for enlarged set of variables
  x_  = (CouNumber *) realloc (x_,  nTotVar * sizeof (CouNumber));
  lb_ = (CouNumber *) realloc (lb_, nTotVar * sizeof (CouNumber));
  ub_ = (CouNumber *) realloc (ub_, nTotVar * sizeof (CouNumber));

  // make expression library point to new vectors
  expression::update (x_, lb_, ub_);

  //print ();

  //for (int i=nVars (), j=0; j < nAuxs (); i++, j++) {
  for (int i=0; i < nVars (); i++) 
    if (variables_ [i] -> Type () == AUX) {

      // initial auxiliary bounds are infinite (they are later changed
      // through branching)

      lb_ [i] = -COUENNE_INFINITY;
      ub_ [i] =  COUENNE_INFINITY;

      // tighten them with propagated bounds
      variables_ [i] -> crossBounds ();

      // and evaluate them
      x_  [i] = (*(variables_ [i] -> Image ())) ();
      lb_ [i] = (*(variables_ [i] -> Lb    ())) ();
      ub_ [i] = (*(variables_ [i] -> Ub    ())) ();

      /*printf (":: %10g [%10g, %10g] [", x_ [i], lb_ [i], ub_ [i]);

      variables_ [i] -> Lb () -> print (); printf (",");
      variables_ [i] -> Ub () -> print (); printf ("]\n");*/
    }

  //graph_ -> print ();
  graph_ -> createOrder ();
  //graph_ -> print ();

  // TODO: fill in numbering structure

  /*for (int k=100; k--;)
    for (int i = 0; i < nVars (); i++)
      if (Var (i) -> Type () == AUX) {
	lb_ [i] = (*(Var (i) -> Lb ())) ();
	ub_ [i] = (*(Var (i) -> Ub ())) ();
	}*/

  //print ();

  delete graph_;
  graph_ = NULL;
}
