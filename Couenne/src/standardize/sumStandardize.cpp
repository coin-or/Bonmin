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

//#define DEBUG


/// given (expression *) element of sum, returns (coe,ind0,ind1)
/// depending on element:
///
/// 1) a * x_i ^ 2   ---> (a,i,?)   return COU_EXPRPOW
/// 2) a * x_i       ---> (a,i,?)   return COU_EXPRVAR
/// 3) a * x_i * x_j ---> (a,i,j)   return COU_EXPRMUL
/// 4) a             ---> (a,?,?)   return COU_EXPRCONST
///
/// x_i and/or x_j may come from standardizing other (linear or
/// quadratic operator) sub-expressions

void decomposeTerm (CouenneProblem *p, expression *term,
		    CouNumber initCoe,
		    CouNumber &c0,
		    std::map <int,                 CouNumber> &lmap,
		    std::map <std::pair <int,int>, CouNumber> &qmap);


/// general procedure to standardize a sum under different forms
/// (exprGroup, exprSum, exprSub, exprOpp)
exprAux *linStandardize (CouenneProblem *, 
			 bool addAux, 
			 CouNumber, 
			 std::map <int,                 CouNumber> &,
			 std::map <std::pair <int,int>, CouNumber> &);

/// translate a sum/difference/exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSum::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  int cod = code ();

  CouNumber c0 = 0; // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // initialize linear/quad maps with the original values/indices of
  // the linear part

  if ((cod == COU_EXPRGROUP) ||
      (cod == COU_EXPRQUAD)) {  // fill linear structure

    exprGroup *eg = dynamic_cast <exprGroup *> (this);

    c0 += eg -> getc0 ();

    int       *olind = eg -> getIndices ();
    CouNumber *olcoe = eg -> getCoeffs ();

    for (int i = eg -> getnLTerms (); i--;)
      lmap.insert (std::pair <int, CouNumber> (olind [i], olcoe [i]));

    if (cod == COU_EXPRQUAD) { // fill quadratic structure

      exprQuad *eq = dynamic_cast <exprQuad *> (this);

      int       *oqindI = eq -> getQIndexI ();
      int       *oqindJ = eq -> getQIndexJ ();
      CouNumber *oqcoe  = eq -> getQCoeffs ();

      for (int i = eq -> getnQTerms (); i--;) {
	std::pair <int, int> ind (oqindI [i], oqindJ [i]);
	qmap.insert (std::pair <std::pair <int, int>, CouNumber> (ind, oqcoe [i]));
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  for (int i=0; i<nargs_; i++)
    decomposeTerm (p, arglist_ [i], 1, c0, lmap, qmap);

  ////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
  printf ("decompTerm returns: [");
  for (std::map <int, CouNumber>::iterator i = lmap.begin (); i != lmap.end (); ++i)
    printf ("<%d,%g>", i -> first, i -> second);
  printf ("] -- [");
  for (std::map <std::pair <int, int>, CouNumber>::iterator i = qmap.begin (); i != qmap.end (); ++i)
    printf ("<%d,%d,%g>", i -> first.first, i -> first.second, i -> second);
  printf ("] (%g)\n", c0);
#endif

  return linStandardize (p, addAux, c0, lmap, qmap);
}


/// translate an exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprOpp::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  CouNumber c0 = 0;   // final constant term

  decomposeTerm (p, argument_, -1, c0, lmap, qmap);

  return linStandardize (p, addAux, c0, lmap, qmap);
}


/// translate a difference (exprSub) into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSub::standardize (CouenneProblem *p, bool addAux) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  CouNumber c0 = 0;   // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  decomposeTerm (p, arglist_ [0],  1, c0, lmap, qmap);
  decomposeTerm (p, arglist_ [1], -1, c0, lmap, qmap);

  return linStandardize (p, addAux, c0, lmap, qmap);
}
