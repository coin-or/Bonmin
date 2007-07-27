/*
 * Name:    analyzeSparsity.cpp
 * Author:  Pietro Belotti
 * Purpose: return one or more exprGroup/exprQuad based on sparsity of
 *          original one
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <map>
#include <CouenneTypes.h>
#include <CouenneProblem.hpp>
#include <exprQuad.hpp>
#include <exprMul.hpp>
#include <exprPow.hpp>

/// analyze sparsity of potential exprQuad/exprGroup and change
/// linear/quadratic maps accordingly, if necessary by adding new
/// auxiliary variables and including them in the linear map
void analyzeSparsity (CouenneProblem *p, CouNumber c0, 
		      std::map <int,                 CouNumber> &lmap,
		      std::map <std::pair <int,int>, CouNumber> &qmap) {

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.begin ();
       i != qmap.end (); i++) {

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
    //if (aux)
    linsert (lmap, aux -> Index (), i -> second);
  }

  qmap.erase (qmap.begin (), qmap.end ());
  /*
  for (std::map <int, CouNumber>::iterator i = lmap.begin ();
       i != lmap.end (); i++) 
    printf ("[%4d %10g]\n", i -> first, i -> second);

  for (std::map <std::pair <int,int>, CouNumber>::iterator i = qmap.begin ();
       i != qmap.end (); i++) 
    printf ("<%4d %4d %10g>\n", i -> first.first, i -> first.second, i -> second);
  */
  return;

  if (qmap.size () == 1) {

    // very simple case: we have a linear term plus a single bilinear
    // x*y (or square x^2) term. 
  }

  // in general, decompose qmap+lmap into several (qmap+lmap)s so that
  // each corresponds to an exprQuad to be transformed into a single
  // auxiliary variable
}
