/*
 * Name:    constrStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardization of constraints
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblemElem.hpp>
#include <CouenneProblem.hpp>


/// replace, in all expressions of the problem (auxiliaries,
/// objectives and constraints) link to an original variable that has
/// gone auxiliary

void CouenneProblem::auxiliarize (exprAux *aux) {

  // find original variable with index = w -> Index ()

  int index = aux -> Index ();

  if (index < 0) {
    printf ("CouenneProblem::auxiliarize (-1)\n");
    return;
  }

  std::vector <exprVar *>::iterator orig;

  for (orig  = variables_.begin ();
       orig != variables_.end (); orig++)

    if ((*orig) -> Index () == index) // found it!
      break;

  if (orig == variables_ . end ()) {
    printf ("CouenneProblem::auxiliarize: no original variable\n");
    return;
  }

  // all objectives

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i)
    (*i) -> Body () -> replace (*orig, aux);

  // and all constraints

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i)
    if ((*i) -> Body ()) 
      (*i) -> Body () -> replace (*orig, aux);

  // substitute it with w in all auxiliaries

  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)
    if (((*i) -> Type () == AUX) && 
	((*i) -> Index () != (*orig) -> Index ()))
      (*i) -> Image () -> replace (*orig, aux);

  // replace it with new auxiliary

  *orig = aux;
}
