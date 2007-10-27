/*
 * Name:    analyzeSparsity.cpp
 * Author:  Pietro Belotti
 * Purpose: return one or more exprGroup/exprQuad based on sparsity of
 *          original one
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <map>
#include <set>

#include "CouenneTypes.hpp"
#include "CouenneProblem.hpp"

#include "exprQuad.hpp"
#include "exprMul.hpp"
#include "exprPow.hpp"

#define MIN_DENSITY 0.5

/// analyze sparsity of potential exprQuad/exprGroup and change
/// linear/quadratic maps accordingly, if necessary by adding new
/// auxiliary variables and including them in the linear map
void analyzeSparsity (CouenneProblem *p, 
		      CouNumber c0, 
		      std::map <int,                 CouNumber> &lmap,
		      std::map <std::pair <int,int>, CouNumber> &qmap) {

  //  return; // comment this if you don't want exprQuad's around

  // simple technique: if number of elements in quadratic map is more
  // than a given fraction of n^2, then turn it into an exprQuad,
  // otherwise break it down. Count n first.

  std::set <int> occur;

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.begin ();
       i != qmap.end (); ++i) {

    if (occur.find (i -> first.first) == occur.end ()) 
      occur.insert (i -> first.first);

    if ((i -> first.first != i -> first.second) &&
	(occur.find (i -> first.second) == occur.end ()))
      occur.insert (i -> first.second);
  }

  /*
  printf ("qmap has %d element, occur has %d, md*s*(s+1)/2 = %g\n", 
	  qmap.size (), 
	  occur.size (),
	  MIN_DENSITY * occur.size () * (occur.size () + 1) / 2);
  */

  if ((qmap.size () > MIN_DENSITY * occur.size () * (occur.size () + 1) / 2) &&
      (occur.size () > 4))

    return;

  // flatten exprQuad to a sum of terms (disaggregate). This is while
  // we are testing exprQuad's

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.begin ();
       i != qmap.end (); ++i) {

    int indI = i -> first.first,
        indJ = i -> first.second;

    exprAux *aux = (indI != indJ) ? 
      p -> addAuxiliary 
      (new exprMul (new exprClone (p -> Var (indI)),
		    new exprClone (p -> Var (indJ)))) : 
      p -> addAuxiliary 
      (new exprPow (new exprClone (p -> Var (indI)),
		    new exprConst (2)));

    //    aux -> print (); printf (" := "); aux -> Image () -> print (); printf ("\n");

    linsert (lmap, aux -> Index (), i -> second);
  }

  if (qmap.size () == 1) {

    // very simple case: we have a linear term plus a single bilinear
    // x*y (or square x^2) term. 
  }

  qmap.erase (qmap.begin (), qmap.end ());

  // in general, decompose qmap+lmap into several (qmap+lmap)s so that
  // each corresponds to an exprQuad to be transformed into a single
  // auxiliary variable

  // build graph and look for components -- TODO
}
