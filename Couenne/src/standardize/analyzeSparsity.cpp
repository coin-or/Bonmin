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
#include "lqelems.hpp"

#define MIN_DENSITY 0.5

//#define DEBUG

/// analyze sparsity of potential exprQuad/exprGroup and change
/// linear/quadratic maps accordingly, if necessary by adding new
/// auxiliary variables and including them in the linear map
void CouenneProblem::analyzeSparsity (CouNumber c0, 
				      LinMap &lmap,
				      QuadMap &qmap) {

  if (qmap.Map().size () == 0) return;

  // simple technique: if number of elements in quadratic map is more
  // than a given fraction of n^2, then turn it into an exprQuad,
  // otherwise break it down. Count n first.

  std::set <int> occur;
  unsigned int nsquares = 0;

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.Map().begin ();
       i != qmap.Map().end (); ++i) {

    int
      first  = i -> first.first,
      second = i -> first.second;

    if (occur.find (first) == occur.end ()) 
      occur.insert (first);

    if (first != second) {
      if  (occur.find (second) == occur.end ())
	occur.insert (second);
    } else nsquares++;
  }

#ifdef DEBUG
  printf ("qmap has %d element, occur has %d, md*s*(s+1)/2 = %g\n", 
	  qmap.Map().size (), 
	  occur.size (),
	  MIN_DENSITY * occur.size () * (occur.size () + 1) / 2);
#endif

  int nterms = occur.size ();

  if ((qmap.Map().size () >= MIN_DENSITY * nterms * (nterms+1) / 2) && (nterms >= 2)
      //|| (nsquares > nterms/2)
      || (nsquares >= occur.size ())
      )
    return; // keep current exprQuad structure

  // flatten exprQuad to a sum of terms (disaggregate). This is while
  // we are testing exprQuad's

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.Map().begin ();
       i != qmap.Map().end (); ++i) {

    int indI = i -> first.first,
        indJ = i -> first.second;

    exprAux *aux = (indI != indJ) ? 
      addAuxiliary 
      (new exprMul (new exprClone (Var (indI)),
		    new exprClone (Var (indJ)))) : 
      addAuxiliary 
      (new exprPow (new exprClone (Var (indI)),
		    new exprConst (2.)));

    //    aux -> print (); printf (" := "); aux -> Image () -> print (); printf ("\n");

    lmap.insert (aux -> Index (), i -> second);
  }

  if (qmap.Map().size () == 1) {

    // very simple case: we have a linear term plus a single bilinear
    // x*y (or square x^2) term. 
  }

  qmap.Map().erase (qmap.Map().begin (), qmap.Map().end ());

  // in general, decompose qmap+lmap into several (qmap+lmap)s so that
  // each corresponds to an exprQuad to be transformed into a single
  // auxiliary variable

  // build graph and look for components -- TODO
}
