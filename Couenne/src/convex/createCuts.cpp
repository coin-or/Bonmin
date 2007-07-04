/*
 * Name:    createCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: a standard cut creator for use with convexification
 *
 * (C) 2006 Pietro Belotti, all rights reserved.
 * This file is distributed under the Common Public License
 */

#include <OsiRowCut.hpp>

#include <CouenneTypes.h>
#include <CouennePrecisions.h>
#include <CouenneCutGenerator.hpp>
#include <CouenneProblem.hpp>


/// general procedure for inserting a linear cut with up to three
/// variables. Return 1 if cut inserted, 0 if none, <0 if error

int CouenneCutGenerator::createCut (OsiCuts &cs,
				    CouNumber rhs, int sign, 
				    int i1, CouNumber c1,
				    int i2, CouNumber c2,
				    int i3, CouNumber c3,
				    bool is_global)       const {
  bool numerics = false;

  static bool warned_large_coeff = false;

  // a maximum of three terms are allowed here. If index -1 (default)
  // the term is not considered

  int nterms = 0;

  if (i1 >= 0) {if (fabs (c1) > COU_MAX_COEFF) numerics = true; nterms++;} else c1 = 0;
  if (i2 >= 0) {if (fabs (c2) > COU_MAX_COEFF) numerics = true; nterms++;} else c2 = 0;
  if (i3 >= 0) {if (fabs (c3) > COU_MAX_COEFF) numerics = true; nterms++;} else c3 = 0;

  if (!nterms) // nonsense cut
    return 0;

  // cut has large coefficients/rhs, bail out
  if (numerics || (fabs (rhs) > COU_MAX_COEFF)) {
    if (!warned_large_coeff) {
      printf ("### Discarding cut, large coeff/rhs: %g (%d), %g (%d), %g (%d); %g\n", 
	      c1, i1, c2, i2, c3, i3, rhs);
      warned_large_coeff = true;
    }
    return 0;
  }

  if (!firstcall_ && addviolated_) { // need to check violation 

    CouNumber *x = const_cast <CouNumber *> (problem_ -> X ());

    // compute violation
    CouNumber violation = - rhs;

    if (i1 >= 0) violation += c1 * x [i1];
    if (i2 >= 0) violation += c2 * x [i2];
    if (i3 >= 0) violation += c3 * x [i3];

    // return 0 if not violated

    if (((violation <   COUENNE_EPS) || (sign > 0)) &&
	((violation > - COUENNE_EPS) || (sign < 0)))
      return 0;
  }

  // check if cut violates optimal solution (irrespective of the
  // branching rules applied, so handle with care)

  CouNumber *best = problem_ -> bestSol (), lhs = 0;

  bool print = false;

  if (best) {

    if (i1 >= 0) lhs += c1 * best [i1];
    if (i2 >= 0) lhs += c2 * best [i2];
    if (i3 >= 0) lhs += c3 * best [i3];

    if ((sign <= 0) && (lhs > rhs + COUENNE_EPS)) 
      {printf ("### cut (%d,%d,%d) (%g,%g,%g) violates optimum: %g >= %g [%g]\n", 
	       i1,i2,i3, c1,c2,c3, lhs, rhs, lhs - rhs); print = true;}

    if ((sign >= 0) && (lhs < rhs - COUENNE_EPS)) 
      {printf ("### cut (%d,%d,%d) (%g,%g,%g) violates optimum: %g <= %g [%g]\n", 
	       i1,i2,i3, c1,c2,c3, lhs, rhs, rhs - lhs); print = true;}
  }

  // You are here if:
  //
  // 1) this is the first call to CouenneCutGenerator::generateCuts()
  // 2) we also want unviolated cuts
  // 3) the cut is violated

  // two cases: cut is of the form w1 [<|>]= alpha, hence a column
  // cut, or it is of the form (a w1 + b w2 + c w3 [<|>]= alpha), a
  // row cut

  if ((i2 < 0) && (i3 < 0)) { // column cut

    if ((fabs (c1) < COUENNE_EPS) && (fabs (rhs) > 1e10*COUENNE_EPS)) {
      printf ("#### nonsense column cut: %e w_%d <>= %e\n", c1, i1, rhs);
      return 0;
    }

    OsiColCut *cut = new OsiColCut;

    CouNumber bound = rhs/c1;

    if (c1 < 0) sign = -sign;

    CouNumber &curL = problem_ -> Lb (i1),
              &curU = problem_ -> Ub (i1);

    if (sign <= 0) {
      cut -> setUbs (1, &i1, &bound); 
      if (bound < curU - COUENNE_EPS) {
	curU = bound;
      }
    }

    if (sign >= 0) {
      cut -> setLbs (1, &i1, &bound); 
      if (bound > curL + COUENNE_EPS) {
	curL = bound;
      }
    }

    cut -> setGloballyValid (is_global); // global?

    //    if (print) cut -> print ();

    cs.insert (cut);
    delete cut;

  } else { // row cut

    CouNumber *coeff = new CouNumber [nterms]; 
    int       *index = new int       [nterms];
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
    //    if (print) cut -> print ();
    cs.insert (cut);
    delete cut;
  }

  return 1;
}
