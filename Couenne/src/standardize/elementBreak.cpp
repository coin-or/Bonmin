/*
 * Name:    elementBreak.cpp
 * Author:  Pietro Belotti
 * Purpose: decompose element of sum if it is of the form cx or x
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblemElem.hpp>
#include <CouenneProblem.hpp>

#include <exprSum.hpp>
#include <exprGroup.hpp>


/// given an element of a sum, check if it is a variable (possibly
/// with a coefficient) and return its index (and the coefficient) if
/// it has not been spotted as an auxiliary (check wentAux)

void elementBreak (expression *arg, int &index, CouNumber &coeff) {

  CouNumber oppMulCoe = 1.;

  bool isMul = false;

  index = -1;

  if (arg -> Linearity () <= LINEAR) { 

    // check if it's w, c*w, or w/c. If more complicated (sum or
    // subtraction) nevermind

    switch (arg -> code ()) { // check element of sum

    case COU_EXPRCONST: break; // it is a constant, nevermind

    case COU_EXPRVAR: // it is a simple variable, w
      index = arg -> Index ();
      coeff = 1.;
      break;

    case COU_EXPROPP: { // it is the opposite of a simple variable or of
                        // another linear term

      index = arg -> Argument () -> Index (); 
      coeff = -1.;

      oppMulCoe = -1; // to be used below
      arg = arg -> Argument (); // transfer analysis to argument of exprOpp

      int argcode = arg -> code ();

      if      (argcode == COU_EXPRMUL) isMul = true;
      else if (argcode != COU_EXPRDIV) break;
    } 
      // no break. Bail out of switch only if this was a -w rather
      // than a -c*w or a -w/c

    case COU_EXPRMUL: 
      if (isMul) { 

      // it is c*w or w*c. Forget the case with more than two
      // non-constant arguments.

      expression **factors = arg -> ArgList ();

      int pos = 0;

      index = (*factors) -> Index ();

      if (index < 0) 
	index = factors [pos = 1] -> Index ();

      if ((index >= 0) && 
	  (fabs (coeff = oppMulCoe * factors [1 - pos] -> Value ()) < COUENNE_EPS))
	index = -1;

      break;
    } // no break outside, there could be some residue from case COU_EXPROPP

    case COU_EXPRDIV: { // if linear, it must be of the form w/c

      expression **factors = arg -> ArgList ();

      coeff = factors [1] -> Value ();
      index = (*factors) -> Index ();

      if (fabs (coeff) < COUENNE_EPS)
	index = -1;
      else coeff = 1. / coeff;

    } break;

    default: break;
    }
  }
}
