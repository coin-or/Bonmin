/*
 * Name:    standardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize all expressions in a problem
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprIVar.hpp"
#include "exprClone.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "depGraph.hpp"

//#define DEBUG

/// standardize (nonlinear) common expressions, objectives, and constraints

void CouenneProblem::standardize () {

  /*printf ("current point: %d vars -------------------\n", domain_.current () -> Dimension ());
  for (int i=0; i<domain_.current () -> Dimension (); i++)
  printf ("%20g [%20g %20g]\n", domain_.x (i), domain_.lb (i), domain_.ub (i));*/

  // create dependence graph to assign an order to the evaluation (and
  // bound propagation, and -- in reverse direction -- implied bounds)
  graph_ = new DepGraph;

  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); ++i)
    graph_ -> insert (*i);

  // allocate space in auxiliaries_ from commonexprs_

  int initVar = variables_ . size () - commonexprs_ . size ();

  // Defined variables ///////////////////////////////////////////

  // standardize initial aux variables (aka defined variables, aka
  // common expression)

#ifdef DEBUG
  if (commonexprs_.size ()) printf ("%d common exprs, initVar = %d = %d - %d\n", 
				    commonexprs_.size (), 
				    initVar, 
				    variables_ . size (), 
				    commonexprs_ . size ());
#endif

  for (std::vector <expression *>::iterator i = commonexprs_ . begin ();
       i != commonexprs_ . end (); ++i) {

#ifdef DEBUG
    printf ("-- stdz common expr [%d] :=", initVar); fflush (stdout);
    (*i) -> print (); printf ("\n"); fflush (stdout);
#endif

    exprAux *naux = (*i) -> standardize (this, false);

    expression *img = naux -> Image ();

    exprAux *newvar = new exprAux (img, initVar, 1 + img -> rank (), exprAux::Unset, &domain_);
    //img -> isInteger () ? exprAux::Integer : exprAux::Continuous);

    auxiliarize (newvar); // takes care of putting newvar at right position in variables_

    graph_ -> insert (newvar);
    graph_ -> erase (naux);

    //variables_ . erase (variables_ . end () - 1);

#ifdef DEBUG
    if (naux) {
      printf ("done: "); fflush (stdout);
      naux -> print (); printf ("\n");
      printf (" := "); fflush (stdout);
      naux -> Image () -> print (); printf ("\n..."); fflush (stdout);
    } else if (*i) {
      (*i) -> print ();
      //printf (" := "); fflush (stdout);
      //aux -> Image () -> print (); 
      printf ("\n");
    } else printf ("[n]aux NULL!\n");
#endif

    initVar++;
  }

  // OBJECTIVES //////////////////////////////////////////////////////////////////////////////

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) {

#ifdef DEBUG
    printf ("Objective "); (*i) -> print ();
#endif

    exprAux *aux = (*i) -> standardize (this);

#ifdef DEBUG
    printf ("      --> "); (*i) -> print ();
#endif

    if (aux) 
      (*i) -> Body (new exprClone (aux));

#ifdef DEBUG
    printf ("      --> "); (*i) -> print (); printf ("...................\n");
#endif
  }

  // commuted_ is an array with a flag for each original variable,
  // which is true at position i if initially original variable x_i
  // went auxiliary

  commuted_ = new bool [nVars ()];
  for (int i = nVars (); i--;)
    *commuted_++ = false;
  commuted_ -= nVars ();

  //  std::vector <CouenneConstraint *> con2;
  std::vector <std::vector <CouenneConstraint *>::iterator> iters2erase;


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
      //      con2.push_back (*i);
    }
    else iters2erase.push_back (i);

    //(*i) -> Body () -> realign (this);

#ifdef DEBUG
    printf (" --> ");
    (*i) -> print ();
    printf ("..............................................................\n");
    //print ();
#endif

    /*printf ("=== "); fflush (stdout); 
    aux -> print (); printf (" := "); fflush (stdout);
    aux -> Image () -> print (); printf ("\n");*/
  }

  for (unsigned int i = iters2erase.size (); i--;)
    constraints_. erase (iters2erase [i]);

  // Look for auxiliaries of the form w:=x and replace each occurrence of w with x

  //if (0)
  for (std::vector <exprVar *>::iterator i = variables_.begin (); 
       i != variables_.end (); ++i)

    if ((*i) -> Type () == AUX) {

      int type = (*i) -> Image () -> Type ();

      if ((type == VAR) || (type == AUX)) {

	// found w_k = x_h. 
	// 
	// Check if either is integer, the survivor will be integer too
	// Replace all occurrences of w_k with x_h

	/*printf ("redundancy: "); 
	(*i)             -> print (); printf (" := "); 
	(*i) -> Image () -> print (); printf ("\n");*/

	// use the info on the variable to be discarded: tighter
	// bounds and integrality that the replacement might not have.

	int 
	  indStays  = (*i) -> Image () -> Index (), // index h
	  indLeaves = (*i)             -> Index (); // index k

	if (indStays == indLeaves)  // very strange case, w_h = x_h
	  continue;

	// do not swap! x_h could be in turn an auxiliary...
	//
	//if (indStays > indLeaves) 
	//{int swap = indStays; indStays = indLeaves; indLeaves = swap;} // swap

	exprVar 
	  *varStays  = variables_ [indStays],
	  *varLeaves = variables_ [indLeaves];

	// intersect features of the two variables (integrality, bounds)

	varStays -> lb () = varLeaves -> lb () = CoinMax (varStays -> lb (), varLeaves -> lb ());
	varStays -> ub () = varLeaves -> ub () = CoinMin (varStays -> ub (), varLeaves -> ub ());

	if (varStays  -> isInteger () ||
	    varLeaves -> isInteger ()) {

	  varStays -> lb () = ceil  (varStays -> lb ());
	  varStays -> ub () = floor (varStays -> ub ());

	  if (varStays -> Type () == AUX)
	    varStays -> setInteger (true);
	  else {
	    //expression *old = varStays; // !!! leak
	    variables_ [indStays] = varStays = new exprIVar (indStays, &domain_);
	    auxiliarize (varStays); // replace it everywhere in the problem
	    //delete old;
	  }
	}

	auxiliarize (varLeaves, varStays); // now replace occurrences of w_k with x_h

	//if (varLeaves -> Index () >= nOrig_) // why check? It's not there anymore.
	varLeaves -> zeroMult ();
      }
    }

  // TODO: re-compute ranks

#ifdef DEBUG
  // Use with caution. Bounds on auxs are not defined yet, so valgrind complains
  printf ("done with standardization:\n");
  //print (); 
#endif

  // Create evaluation order ////////////////////////////////////////////////////

  delete auxSet_;

  // reallocate space for enlarged set of variables
  domain_.current () -> resize (nVars ());

  //graph_ -> print ();
  graph_ -> createOrder ();
  //graph_ -> print ();

  assert (graph_ -> checkCycles () == false);

  // Fill numbering structure /////////////////////////////////////////////////

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

      if (variables_ [ord] -> Index () >= nOrig_) { // and one that was not an original, originally...

	domain_.lb (ord) = -COIN_DBL_MAX;
	domain_.ub (ord) =  COIN_DBL_MAX;
      }

      //printf ("from "); variables_ [ord] -> Lb    () -> print (); 

      // tighten them with propagated bounds
      variables_ [ord] -> crossBounds ();

      //printf ("to "); variables_ [ord] -> Lb    () -> print (); printf (", now eval\n");

      // and evaluate them
      domain_.x  (ord) = (*(variables_ [ord] -> Image ())) ();
      domain_.lb (ord) = (*(variables_ [ord] -> Lb    ())) ();
      domain_.ub (ord) = (*(variables_ [ord] -> Ub    ())) ();

#ifdef DEBUG
      printf (":::: %10g [%10g, %10g] [", domain_.x (ord), domain_.lb (ord), domain_.ub (ord));
      variables_ [ord] -> Lb () -> print (); printf (",");
      variables_ [ord] -> Ub () -> print (); printf ("]\n");
#endif

      bool integer = variables_ [ord] -> isInteger ();

      if (integer) {
	domain_.lb (ord) = ceil  (domain_.lb (ord) - COUENNE_EPS);
	domain_.ub (ord) = floor (domain_.ub (ord) + COUENNE_EPS);
      }
    }
  }

  delete [] commuted_;  commuted_ = NULL;
  delete    graph_;     graph_    = NULL;
}
