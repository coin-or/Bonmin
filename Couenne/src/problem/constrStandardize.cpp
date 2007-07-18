/*
 * Name:    constrStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization of constraints
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblemElem.hpp>
#include <CouenneProblem.hpp>

#include <exprAux.hpp>
#include <depGraph.hpp>


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

  //printf ("||||||| standardizing constraint: "); print ();

  /*printf (" ["); fflush (stdout);
  lb_ -> print ();
  printf (","); fflush (stdout);
  ub_ -> print ();
  printf ("] \n");*/

  if (0) // Back to normal
  if (compareExpr (&lb_, &ub_) == 0) { // this is an equality constraint

    expression *rest;

    // split w from f(x)
    int wind = splitAux (p, (*lb_) (), body_, rest, p -> Commuted ());

    if (wind >= 0) { // this IS the definition of an auxiliary variable w = f(x)

      p -> Commuted () [wind] = true;

      //printf ("---> %d & ", wind); fflush (stdout);
      //rest -> print (); printf ("\n");

      exprAux *w = new exprAux (rest, wind, 1 + rest -> rank (p));

      // no such expression found in the set:
      p -> Auxiliaries () . push_back (w); // 1) create entry therein
      p -> AuxSet      () -> insert   (w); // 2) beware of useless copies
      p -> getDepGraph () -> insert   (w); // 3) introduce it in acyclic structure

      // replace ALL occurrences of original variable (with index
      // wind) with newly created auxiliary
      p -> auxiliarize (w);

      return NULL;
    }
  }

  //printf ("\nnormal\n-----------------\n");

  return body_ -> standardize (p);
}
