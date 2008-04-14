/*
 * Name:    auxiliarize.cpp
 * Author:  Pietro Belotti
 * Purpose: replace occurrences of original variable in a problem with
 *          auxiliary with the same index
 *
 * (C) Carnegie-Mellon University, 2007-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "exprClone.hpp"

/// replace, in all expressions of the problem (auxiliaries,
/// objectives and constraints) link to an original variable that has
/// gone auxiliary

void CouenneProblem::auxiliarize (exprVar *aux, exprVar *subst) {

  bool same_var = (subst == NULL);

  if (!subst) 
    subst = aux;

  // find original variable with index = w -> Index ()

  int index = aux -> Index ();

  assert (index >= 0);

  std::vector <exprVar *>::iterator orig;

  for (orig  = variables_.begin ();
       orig != variables_.end (); ++orig)

    if ((((*orig) -> Type () == VAR) || !same_var) && 
	((*orig) -> Index () == index)) // found it

      break;

  if (orig == variables_ . end ()) {
    printf ("CouenneProblem::auxiliarize: no original variables correspond\n");
    return;
  }

  // all objectives

  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) {

    expression *body = (*i) -> Body ();

    if (body) {
      if ((body -> Type () == VAR) || 
	  (body -> Type () == AUX)) {

	if (body -> Index () == (*orig) -> Index ()) {
      
	  delete body;//(*i) -> Body ();
	  (*i) -> Body (new exprClone (subst));
	}
      } else body -> replace (*orig, subst);
    }
  }

  // and all constraints

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i) {

    expression *body = (*i) -> Body ();

    if (body) {
      if ((body -> Type () == VAR) ||
	  (body -> Type () == AUX)) {

	if (body -> Index () == (*orig) -> Index ()) {
      
	  delete body;//(*i) -> Body ();
	  (*i) -> Body (new exprClone (subst));
	}
      } else body -> replace (*orig, subst);
    }
  }

  // substitute it with w in all auxiliaries

  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)

    if (((*i) -> Type () == AUX) &&                  // replace in all aux's image
	((*i) -> Index () != (*orig) -> Index ())) { // skip same variable

      expression *image = (*i) -> Image ();

      if ((image -> Type () == VAR) || 
	  (image -> Type () == AUX)) {
	if (image -> Index () == (*orig) -> Index ()) {

	  //delete image;//(*i) -> Image (); // !!! Possible leak? TODO
	  (*i) -> Image (new exprClone (subst));
	} //else (*i) -> Image () -> replace (*orig, aux);
      } else image  -> replace (*orig, subst);
    }

  // replace it with new auxiliary

  if (same_var)
    *orig = aux;
}
