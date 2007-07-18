/*
 * Name:    quadCuts.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: based on upper and lower convexification, add cuts to convexify
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>


/// [Pierre] Given the data in dCoeff_ and dIndex_, add
/// convexification cuts to the OsiCuts structure
void exprQuad::quadCuts (OsiCuts &cs, CouenneCutGenerator *cg) {

  /// retrieve linear coefficients/indices and constant from member of
  /// the base class, exprGroup

  int  nterms = nlterms_ + nargs_,
      *index  = new int [nterms];

  CouNumber *coeff   = new CouNumber [nterms],
             a0      = c0_;

  // linear part of exprGroup
  for (nterms = 0; nterms < nlterms_; nterms++) {
    index [nterms] = index_ [nterms];
    coeff [nterms] = coeff_ [nterms];
  }

  // deal with nonlinear part (has been standardized into sum of
  // auxiliary variables)

  for (int i=0; i<nargs_; i++) {

    int ind = arglist_ [i] -> Index ();

    if (ind >= 0) {

      index [nterms]   = ind;
      coeff [nterms++] = 1.;

    } else 
      if (arglist_ [i] -> Type () == CONST) // if not a constant, quit
	a0 += arglist_ [i] -> Value ();
      else {
	printf ("non constant, non variable term in standardized exprSum\nAborting.\n");
	exit (-1);
      }
  }

  // Pierre: now index[], coeff[], and a0 contain the linear data of
  // expression a0 + a^Tx + x^T Q x. Variable nterms contains the
  // number of linear terms.
  //
  // The convexification cuts should take into account the linear data
  // as well, so if we have an auxiliary variable
  //
  // w = a0 + ax + x'Qx
  //
  // the relative lower envelope cut should be of the form 
  //
  // w >= a0 + ax + b + bx
  //
  // where w >= b + bx is the convexification cut of the quadratic
  // expression alone

  delete [] coeff;
  delete [] index;
}
