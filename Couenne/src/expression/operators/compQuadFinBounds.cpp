/*
 * Name:    compQuadFinBounds.cpp
 * Author:  Pietro Belotti
 * Purpose: compute bounds on quadratic form without the contribution of a single variable
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>


void updateInf (CouNumber, CouNumber, CouNumber, int &, int &, int );


void computeQuadFiniteBound (const exprQuad *e,
			     CouNumber &qMin, CouNumber &qMax, 
			     CouNumber *l, CouNumber *u,
			     int &indInfLo, int &indInfUp) {

  exprQuad *eq = const_cast <exprQuad *> (e);

  int 
    nlt = eq -> getnLTerms (),
    *li = eq -> getIndices (),

    nqt = eq -> getnQTerms (),
    *qi = eq -> getQIndexI (),
    *qj = eq -> getQIndexJ ();

  CouNumber
    *lc = eq -> getCoeffs  (),
    *qc = eq -> getQCoeffs ();

  // linear terms ///////////////////////////////////////////////////

  while (nlt--) {

    int ind = *li++;

    CouNumber 
      coe = *lc++, 
      li  = l [ind], 
      ui  = u [ind];

    if (coe < 0) { // negative coefficient

      if (indInfLo > -2) {
	if (ui >  COUENNE_INFINITY) indInfLo = (indInfLo == -1) ? ind : -2;
	else qMin += coe * ui;
      }

      if (indInfUp > -2) {
	if (li < -COUENNE_INFINITY) indInfUp = (indInfUp == -1) ? ind : -2;
	else qMax += coe * li;
      }
    } else { // positive coefficient

      if (indInfLo > -2) {
	if (li < -COUENNE_INFINITY) indInfLo = (indInfLo == -1) ? ind : -2;
	else qMin += coe * li;
      }

      if (indInfUp > -2) {
	if (ui >  COUENNE_INFINITY) indInfUp = (indInfUp == -1) ? ind : -2;
	else qMax += coe * ui;
      }
    }
  }

  // quadratic terms ////////////////////////////////////////////////

  while (nqt--) {

    int i = *qi++,
        j = *qj++;

    CouNumber 
      coe = *qc++,
      lbi = l [i],
      ubi = u [i];

    if (i==j) { // term of the form q_{ii} x_i^2

      CouNumber
	tmin = (ubi <= 0) ? (ubi * ubi) : (lbi >= 0) ? (lbi * lbi) : 0,
	tmax = mymax (lbi*lbi, ubi*ubi);

      if (coe < 0) { // negative coefficient, q_{ii} < 0
	qMax += coe * tmin;
	if (indInfLo > -2) {
	  if (tmax > COUENNE_INFINITY) indInfLo = (indInfLo == -1) ? i : -2;
	  else qMin += coe * tmax;
	}
      } else { // positive coefficient
	qMin += coe * tmin;
	if (indInfUp > -2) {
	  if (tmax > COUENNE_INFINITY) indInfUp = (indInfUp == -1) ? i : -2;
	  else qMax += coe * tmax;
	}
      }
    } else { // term of the form q_{ij} x_i x_j, j\neq i

      CouNumber
	lbj = l [j],
	ubj = u [j];

      coe *= 2;

      // check if index i wrings unbounded value in both directions

      if (lbi < -COUENNE_INFINITY) updateInf (coe, lbj, ubj, indInfUp, indInfLo, i);
      if (lbj < -COUENNE_INFINITY) updateInf (coe, lbi, ubi, indInfUp, indInfLo, j);

      if (ubi >  COUENNE_INFINITY) updateInf (coe, lbj, ubj, indInfLo, indInfUp, i);
      if (ubj >  COUENNE_INFINITY) updateInf (coe, lbi, ubi, indInfLo, indInfUp, j);

      CouNumber term, 
	b1 = coe * lbi * lbj,
	b2 = coe * lbi * ubj,
	b3 = coe * ubi * lbj,
	b4 = coe * ubi * ubj;

      if ((term = mymin (mymin (b1, b2), mymin (b3, b4))) > -COUENNE_INFINITY) qMin += term;
      if ((term = mymax (mymax (b1, b2), mymax (b3, b4))) <  COUENNE_INFINITY) qMax += term;
    }
  }
}


void updateInf (CouNumber coe, CouNumber lb, CouNumber ub, int &indInf1, int &indInf2, int i) {

  if (coe > 0) {
    if (lb < 0) indInf1 = (indInf1 == -1) ? i : -2;
    if (ub > 0) indInf2 = (indInf2 == -1) ? i : -2;
  } else {
    if (lb < 0) indInf2 = (indInf2 == -1) ? i : -2;
    if (ub > 0) indInf1 = (indInf1 == -1) ? i : -2;
  }
}
