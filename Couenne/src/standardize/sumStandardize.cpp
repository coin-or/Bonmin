/*
 * Name:    sumStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: check if expr{Group,Sum,Sub} contains a lot of quadratic/bilinear terms
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprSum.hpp"
#include "exprSub.hpp"
#include "exprOpp.hpp"
#include "exprGroup.hpp"
#include "exprQuad.hpp"
#include "lqelems.hpp"

//#define DEBUG


/// translate a sum/difference/exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSum::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  LinMap lmap;
  QuadMap qmap;

  int cod = code ();

  CouNumber c0 = 0; // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // initialize linear/quad maps with the original values/indices of
  // the linear part

  if ((cod == COU_EXPRGROUP) ||
      (cod == COU_EXPRQUAD)) {  // fill linear structure

    exprGroup *eg = dynamic_cast <exprGroup *> (this);
    exprGroup::lincoeff &lcoe = eg -> lcoeff ();

    c0 += eg -> getc0 ();

    for (exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el)
      lmap.insert (el -> first -> Index (), el -> second);

    if (cod == COU_EXPRQUAD) { // fill quadratic structure

      exprQuad *eq = dynamic_cast <exprQuad *> (this);
      exprQuad::sparseQ &M = eq -> getQ ();

      // derive quadratic part (obtain linear part)
      for (exprQuad::sparseQ::iterator row = M.begin (); row != M.end (); ++row) {

	int xind = row -> first -> Index ();

	for (exprQuad::sparseQcol::iterator col = row -> second.begin (); 
	     col != row -> second.end (); ++col)
	  qmap.insert (xind, col -> first -> Index (), col -> second);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  for (int i=0; i<nargs_; i++)
    p -> decomposeTerm (arglist_ [i], 1, c0, lmap, qmap);

  ////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
  printf ("decompTerm returns: [");
  for (std::map <int, CouNumber>::iterator i = lmap.Map().begin (); i != lmap.Map().end (); ++i)
    printf ("<%d,%g>", i -> first, i -> second);
  printf ("] -- [");
  for (std::map <std::pair <int, int>, CouNumber>::iterator i = qmap.Map ().begin (); 
       i != qmap.Map ().end (); ++i)
    printf ("<%d,%d,%g>", i -> first.first, i -> first.second, i -> second);
  printf ("] (%g)\n", c0);
#endif

  return p -> linStandardize (addAux, c0, lmap, qmap);
}


/// translate an exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprOpp::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  LinMap lmap;
  QuadMap qmap;

  CouNumber c0 = 0;   // final constant term

  p -> decomposeTerm (argument_, -1., c0, lmap, qmap);

  return p -> linStandardize (addAux, c0, lmap, qmap);
}


/// translate a difference (exprSub) into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSub::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  LinMap lmap;
  QuadMap qmap;

  CouNumber c0 = 0;   // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  p -> decomposeTerm (arglist_ [0],  1, c0, lmap, qmap);
  p -> decomposeTerm (arglist_ [1], -1, c0, lmap, qmap);

  return p -> linStandardize (addAux, c0, lmap, qmap);
}
