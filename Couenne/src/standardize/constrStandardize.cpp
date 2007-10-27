/*
 * Name:    constrStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization of constraints
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblemElem.hpp"
#include "CouenneProblem.hpp"

#include "exprAux.hpp"
#include "depGraph.hpp"

//#define DEBUG

/// split a constraint w - f(x) = c into w's index (it is returned)
/// and rest = f(x) + c

int splitAux (CouenneProblem *, CouNumber, expression *, expression *&, bool *);


/// decompose body of constraint through auxiliary variables

exprAux *CouenneConstraint::standardize (CouenneProblem *p) {

  // spot an auxiliary variable in constraint's body w - f(x) and move
  // the explicit w into the vector of auxiliary variables
  //
  // do it only if this is an equality constraint and there is at
  // least one variable that did not show up so far (need a problem
  // structure)

#ifdef DEBUG
  printf ("################################\nStandardizing constraint: "); print ();

  printf (" ["); fflush (stdout);
  lb_ -> print ();
  printf (","); fflush (stdout);
  ub_ -> print ();
  /*  printf ("] {with auxset = ");
  for (std::set <exprAux *, compExpr>::iterator i = p -> AuxSet () -> begin ();
       i != p -> AuxSet () -> end (); i++) {
    printf ("<"); (*i) -> print (); 
    printf (","); (*i) -> Image () -> print (); printf ("> ");
    }*/

  printf ("}\n");
#endif

  if (compareExpr (&lb_, &ub_) == 0) {

    // this is an equality constraint, and as such it could be the
    // definition of an auxiliary

    expression *rest;

    // split w from f(x)
    int wind = splitAux (p, (*lb_) (), body_, rest, p -> Commuted ());

    if (wind >= 0) { // this IS the definition of an auxiliary variable w = f(x)

      p -> Commuted () [wind] = true;

#ifdef DEBUG
      printf ("---> %d & ", wind); fflush (stdout);
      rest -> print (); printf ("\n");
#endif

      // create new variable, it has to be integer if original variable was integer
      exprAux *w = new exprAux (rest, wind, 1 + rest -> rank (p), 
				p -> Var (wind) -> isInteger ());

      std::set <exprAux *, compExpr>::iterator i = p -> AuxSet () -> find (w);

      // no such expression found in the set:
      if (i == p -> AuxSet () -> end ()) {

	p -> AuxSet      () -> insert (w); // 1) beware of useless copies
	p -> getDepGraph () -> insert (w); // 2) introduce it in acyclic structure

#ifdef DEBUG
	printf ("now replacing x [%d] with ", wind);
	w -> print (); printf (" := ");
	w -> Image () -> print (); printf ("\n");
#endif

	// replace ALL occurrences of original variable (with index
	// wind) with newly created auxiliary
	p -> auxiliarize (w);
      } 

#ifdef DEBUG
      else {
	printf ("found aux occurrence of "); 
	w -> print (); printf (" := ");
	w -> Image () -> print (); printf (" ... ");
	(*i) -> print (); printf (" := ");
	(*i) -> Image () -> print (); printf ("\n");
      }
#endif

      return NULL;
    }
  }

#ifdef DEBUG
  printf ("\nnormal\n-----------------\n");
#endif

  return body_ -> standardize (p);
}
