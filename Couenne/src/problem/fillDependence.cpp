/*
 * Name:    fillDependence.cpp
 * Author:  Pietro Belotti
 * Purpose: fill in inverse dependence structure, for CouenneObject
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include <set>

#include "CouenneProblem.hpp"


/// fill in inverse dependence structure: for each variable x give set
/// of auxiliary variables (or better, their indices) whose expression
/// depends on x

void CouenneProblem::fillDependence (Bonmin::BabSetupBase *base) {

  // initialize vector of empty sets
  for (int i=nVars (); i--;)
    dependence_.push_back (std::set <int> ());

  // empty object to fill space for linear-defined auxiliaries and for
  // originals
  CouenneObject nullObject;

  for (std::vector <exprVar *>::iterator i = variables_.begin (); 
       i != variables_.end (); ++i) {

    if (((*i) -> Type () == AUX)                           // consider aux's only
	&& ((*i) -> Image () -> Linearity () > LINEAR)) {  // and nonlinear

      CouenneObject *infeasObj = (*i) -> properObject (this, base, jnlst_);

      if (!infeasObj) // found something that will never be infeasibl
	continue;

      CouenneObject infObj (*infeasObj);

      delete infeasObj;

      // add object for this variable
      objects_.push_back (infObj);

      std::set <int> deplist;

      // fill the set of independent variables on which the expression
      // associated with *i depends; if empty (should not happen...), skip
      if ((*i) -> Image () -> DepList (deplist, STOP_AT_AUX) == 0)
	continue;

      // build dependence set for this variable
      for (std::set <int>::iterator j = deplist.begin (); j != deplist.end (); ++j) {

	std::set <int> &obj = dependence_ [*j];
	int ind = (*i) -> Index ();
	if (obj.find (ind) == obj.end ())
	  obj.insert (ind);
      }

    } else objects_.push_back (nullObject); // null object for original and linear auxiliaries
  }
}
