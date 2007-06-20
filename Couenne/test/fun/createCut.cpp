/*
 * Name:    createCut.cpp
 * Author:  Pietro Belotti
 * Purpose: create simple OsiRowCut (used in initial formulation)
 *
 * (C) Pietro Belotti. This file is published under the Common Public License.
 */

#include <OsiCuts.hpp>

#define COUENNE_EPS 1e-7

/// create a cut. return 1 if cut was inserted into cs, 0 otherwise

int createCut (OsiCuts &cs,           // set of cuts
	       double rhs,            // rhs 
	       int sign,              // -1: $\le$, +1: $\ge$, 0: =
	                              // indices, coefficients 
	                              // (index == -1 means don't care)
	       int i1, double c1,     // first 
	       int i2, double c2,     // second
	       int i3, double c3,     // third
	       bool is_global) {      // is this cut global (true) or local (false)?

  int nterms = 0;

  if (i1 >= 0) nterms++; else c1 = 0;
  if (i2 >= 0) nterms++; else c2 = 0;
  if (i3 >= 0) nterms++; else c3 = 0;

  if (!nterms) // nonsense cut
    return 0;

  double    *coeff = new double [nterms]; 
  int       *index = new int    [nterms];
  OsiRowCut *cut   = new OsiRowCut;

  if (i1 >= 0) {coeff [0] = c1; index [0] = i1;}
  if (i2 >= 0) {coeff [1] = c2; index [1] = i2;}
  if (i3 >= 0) {coeff [2] = c3; index [2] = i3;}

  if (sign <= 0) cut -> setUb (rhs);
  if (sign >= 0) cut -> setLb (rhs);

  cut -> setRow (nterms, index, coeff);

  delete [] coeff;
  delete [] index;

  cut -> setGloballyValid (is_global); // global?
  cs.insert (cut);
  delete cut;

  return 1;
}
