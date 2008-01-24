/*
 * Name:    branchExprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for quadratic forms
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprQuad.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

//#define DEBUG

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprQuad::selectBranch (const CouenneObject *obj, 
				  const OsiBranchingInformation *info,
				  expression *&var, 
				  double * &brpts, 
				  int &way) {

  int ind = -1;

  // use a combination of eigenvectors and bounds

  CouNumber delta = (*(obj -> Reference ())) () - (*this) ();

  /*printf ("infeasibility: ");
  obj -> Reference () -> print (); 
  printf (" [%g=%g] := ", 
	  (*(obj -> Reference ())) (), info -> solution_ [obj -> Reference () -> Index ()]);
	  print (); printf (" [%g]\n", (*this) ());*/

  brpts = (double *) realloc (brpts, sizeof (double));
  way = TWO_RAND;

  // depending on where the current point is w.r.t. the curve,
  // branching is done on the eigenvector corresponding to the minimum
  // or maximum eigenvalue

  std::vector <std::pair <CouNumber, 
    std::vector <std::pair <exprVar *, CouNumber> > > >::iterator         fi = eigen_.begin ();

  std::vector <std::pair <CouNumber, 
    std::vector <std::pair <exprVar *, CouNumber> > > >::reverse_iterator ri = eigen_.rbegin ();

  CouNumber max_span = -COUENNE_INFINITY;

  bool changed_sign = false;

  /////////////////////////////////////////////////////////

  for (;((delta < 0.) && (fi != eigen_. end  ()) || // && (fi -> first < 0.) ||
	 (delta > 0.) && (ri != eigen_. rend ()));  // && (ri -> first > 0.));
       ++fi, ++ri) {

    std::vector <std::pair <exprVar *, CouNumber> > &ev = 
      (delta < 0.) ? fi -> second : ri -> second;

    if ((delta < 0.) && (fi -> first > 0.) ||
	(delta > 0.) && (ri -> first < 0.)) {

      if (max_span > 0.) break; // if found a variable already, return
      changed_sign = true;      // otherwise, keep in mind we are on
				// the wrong eigenvalues' sign
    }

    for (std::vector <std::pair <exprVar *, CouNumber> >::iterator j = ev.begin ();
	 j != ev.end (); ++j) {

      int index = j -> first -> Index ();

      CouNumber 
	lb = info -> lower_ [index],
	ub = info -> upper_ [index];

      // only accept variable if not fixed
      if (fabs (ub-lb) > COUENNE_EPS) {

	  // no variable was found but the eigenvalue is already
	  // positive (negative)
	  //	  changed_sign &&
	  //	  (max_span < 0.))

	CouNumber span = -1;

	if ((lb < -COUENNE_INFINITY) ||
	    (ub >  COUENNE_INFINITY) ||
	    ((span = (ub-lb) * fabs (j -> second)) > max_span + COUENNE_EPS)) {

	  ind = index;
	  var = j -> first;

	  *brpts = midInterval (info -> solution_ [index], lb, ub);

	  if (changed_sign) 
	    break;

	  if (span >= 0) max_span = span; // finite, keep searching
	  else break;                     // span unbounded, stop searching
	}
      }
    }
  }

  if ((eigen_.size () == 0) // if no eigenvalue has been computed yet
      || (ind < 0)) {       // or no index was found, pick largest interval

    CouNumber max_span = -COUENNE_INFINITY;

    for (std::map <exprVar *, std::pair <CouNumber, CouNumber> >::iterator i = bounds_. begin ();
	 i != bounds_. end (); ++i) {

      CouNumber
	lb = i -> second.first,
	ub = i -> second.second,
	span = ub - lb;

      if ((span > COUENNE_EPS) && 
	  (span > max_span)) {

	var = i -> first;
	ind = var -> Index ();
      }
    }

    if (ind < 0) {

      var = obj -> Reference ();
      ind = var -> Index ();

      *brpts = midInterval (info -> solution_ [ind],
			    info -> lower_ [ind],
			    info -> upper_ [ind]);
    }
    else *brpts = midInterval (info -> solution_ [ind], 
			       info -> lower_ [ind],
			       info -> upper_ [ind]);	  

    return fabs (delta);
  }

  return fabs (delta);
}
