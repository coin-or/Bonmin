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


//#define DEBUG

/// standardize (nonlinear) common expressions, objectives, and constraints

void CouenneProblem::standardize () {

  // create dependence graph to assign an order to the evaluation (and
  // bound propagation, and -- in reverse direction -- implied bounds)
  graph_ = new DepGraph;

  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); ++i)
    graph_ -> insert (*i);

  // allocate space in auxiliaries_ from commonexprs_

  int initVar = variables_ . size () - commonexprs_ . size ();

  // DEFINED VARIABLES ///////////////////////////////////////////////////////////////////////

  // standardize initial aux variables (aka defined variables, aka
  // common expression)

  for (int nc = commonexprs_ . size (), i=0; i<nc; i++) {

    expression *aux = commonexprs_ [i];

#ifdef DEBUG
    printf ("////////////// now attempting to standardize defVar "); fflush (stdout);
    aux -> print ();
    //    printf (" := "); fflush (stdout);
    //    aux -> print (); 
    printf ("\n ----> "); fflush (stdout);
#endif

    exprAux *naux = aux -> standardize (this);

    expression *img = naux -> Image ();

    exprAux *newvar = new exprAux (img, initVar, 1 + img -> rank (this));
    variables_ [initVar] = newvar;

    graph_ -> insert (newvar);

    initVar++;

    graph_ -> erase (naux);

    variables_ . erase (variables_ . end () - 1);

#ifdef DEBUG
    if (naux) {
      printf ("done: "); fflush (stdout);
      naux -> print (); printf ("\n");
      printf (" := "); fflush (stdout);
      naux -> Image () -> print (); printf ("\n..."); fflush (stdout);
    } else if (aux) {
      aux -> print ();
      //printf (" := "); fflush (stdout);
      //aux -> Image () -> print (); 
      printf ("\n");
    } else printf ("[n]aux NULL!\n");
#endif
  }

  // OBJECTIVES //////////////////////////////////////////////////////////////////////////////

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) {

#ifdef DEBUG
    printf ("Objective ");
    (*i) -> print ();
#endif

    exprAux *aux = (*i) -> standardize (this);

#ifdef DEBUG
    printf ("      --> ");
      (*i) -> print ();
#endif

    if (aux) 
      (*i) -> Body (new exprClone (aux));

#ifdef DEBUG
    printf ("      --> ");
    (*i) -> print ();
    printf ("...................\n");
#endif
  }

  // commuted_ is an array with a flag for each original variable,
  // which is true at position i if initially original variable x_i
  // went auxiliary

  commuted_ = new bool [nVars ()];
  for (int i = nVars (); i--;)
    *commuted_++ = false;
  commuted_ -= nVars ();

  std::vector <CouenneConstraint *> con2;

  // CONSTRAINTS /////////////////////////////////////////////////////////////////////////////

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin (); 
       i != constraints_.end (); ++i) {

#ifdef DEBUG
    printf ("############# Constraint ");
      (*i) -> print ();
#endif

    exprAux *aux = (*i) -> standardize (this);

#ifdef DEBUG
    printf (" ==> [%d] ", aux ? (aux -> Index ()) : -1);
      (*i) -> print ();
#endif

    if (aux) { // save if standardized
      (*i) -> Body (new exprClone (aux));
      con2.push_back (*i);
    }

#ifdef DEBUG
    printf (" --> ");
    (*i) -> print ();
    printf ("...................\n");
#endif

    /*printf ("=== "); fflush (stdout); 
    aux -> print (); printf (" := "); fflush (stdout);
    aux -> Image () -> print (); printf ("\n");*/
  }

  //  constraints_. erase (constraints_.begin (), constraints_.end ());

  constraints_ = con2;

#ifdef DEBUG
  printf ("done with standardization:\n");
  print ();
#endif

  // CREATE EVALUATION ORDER /////////////////////////////////////////////////////////////////

  delete auxSet_;

  int nTotVar = nVars ();

  // reallocate space for enlarged set of variables
  x_  = (CouNumber *) realloc (x_,  nTotVar * sizeof (CouNumber));
  lb_ = (CouNumber *) realloc (lb_, nTotVar * sizeof (CouNumber));
  ub_ = (CouNumber *) realloc (ub_, nTotVar * sizeof (CouNumber));

  // make expression library point to new vectors
  expression::update (x_, lb_, ub_);

  //  graph_ -> print ();
  graph_ -> createOrder ();
  //  graph_ -> print ();

  // fill numbering structure /////////////////////////////////////////////////

  int n = nVars ();
  numbering_ = new int [n];
  std::set <DepNode *, compNode> vertices = graph_ -> Vertices ();

  for (std::set <DepNode *, compNode>::iterator i = vertices.begin ();
       i != vertices.end (); ++i)
    numbering_ [(*i) -> Order ()] = (*i) -> Index (); 

  //////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < nVars (); i++) {

    int ord = numbering_ [i];

    if (variables_ [ord] -> Type () == AUX) {

      // initial auxiliary bounds are infinite (they are later changed
      // through branching)

      lb_ [ord] = -COUENNE_INFINITY;
      ub_ [ord] =  COUENNE_INFINITY;

      // tighten them with propagated bounds
      variables_ [ord] -> crossBounds ();

      // and evaluate them
      x_  [ord] = (*(variables_ [ord] -> Image ())) ();
      lb_ [ord] = (*(variables_ [ord] -> Lb    ())) ();
      ub_ [ord] = (*(variables_ [ord] -> Ub    ())) ();

#ifdef DEBUG
      printf (":::: %10g [%10g, %10g] [", x_ [ord], lb_ [ord], ub_ [ord]);

      variables_ [ord] -> Lb () -> print (); printf (",");
      variables_ [ord] -> Ub () -> print (); printf ("]\n");
#endif
    }
  }

  //for (int i=0; i<n; i++)
  //printf ("[%4d %4d]\n", i, numbering_ [i]);

  delete graph_;
  graph_ = NULL;
}
