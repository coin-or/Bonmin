/*
 * Name:    conv-exprPow-getBounds.C
 * Author:  Pietro Belotti
 * Purpose: method to get lower and upper bounds of a power x^y
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <CouenneTypes.h>
#include <exprPow.h>
#include <exprConst.h>
#include <exprCopy.h>
#include <exprClone.h>
#include <exprMax.h>
#include <exprMin.h>
#include <exprOpp.h>
#include <CouennePrecisions.h>
#include <CouenneProblem.h>


// compute expressions for lower and upper bounds of a power x^y,
// based on the lower/upper bounds of x and y

void exprPow::getBounds (expression *&lb, expression *&ub) {

  // We have a standardized expression of the form w = x^y, where x or
  // y could be constant. Let us study each case separately.

  if (arglist_ [0] -> Type () == CONST) {

    ////////////////////////////////////////////////////////////////
    // Should already be dealt with by simplify () / standardize ()
    //

    fprintf (stderr, "exprPow::getBounds(): warning, expression is of the type k^x\n");
    return;

  }       /********** end [if (x is constant)] ************/
  else {

    // x is not constant, so it has (possibly different) lower and
    // upper bounds. The expression is x^b, b constant (the case x^y
    // has been decomposed by simplify() into exp(y log x).

    expression *lbbase, *ubbase;
    arglist_ [0] -> getBounds (lbbase, ubbase);

    if (arglist_ [1] -> Type () == CONST) { 

      // expression = x^b, b!=0. There are four cases:
      // 
      // 1) b   is integer and odd  (cube, x^5, etc)
      // 2) b   is integer and even (square, x^8, etc)
      // 3) 1/b is integer and odd  (cube root, x^(1/7), etc)
      // 4) 1/b is integer and even (square root, x^(1/4), etc)
      // 5) none of the above
      //
      // For all of these, need to check if the exponent is negative...

      CouNumber expon = arglist_ [1] -> Value ();
      int rndexp;

      bool isInt    =  fabs (expon - (rndexp = FELINE_round (expon))) < COUENNE_EPS;
      bool isInvInt = !isInt &&  
                      ((fabs (expon) > COUENNE_EPS) && 
		       (fabs (1/expon - (rndexp = FELINE_round (1/expon))) < COUENNE_EPS));

      if ((isInt || isInvInt) && (rndexp % 2) && (rndexp > 0)) { 

	// the exponent is integer (or inverse integer), odd and
	// positive, hence the function is monotone non decreasing

	lb = new exprPow (lbbase, new exprConst (expon));
	ub = new exprPow (ubbase, new exprConst (expon));
      } 
      else {

	// the exponent is either negative, integer even, or fractional

	expression **all = new expression * [6];

	all [0] = new exprOpp   (lbbase);
	all [2] = new exprConst (0);
	all [4] = ubbase;

	if (expon > COUENNE_EPS) 
	     all [1] = new exprPow (new exprCopy (lbbase), new exprConst (expon));
	else all [1] = new exprPow (new exprCopy (ubbase), new exprConst (expon));

	// all [3] is lower bound when lbbase <= 0 <= ubbase

	if (expon > COUENNE_EPS) all [3] = new exprConst (0);
	else if (isInt || isInvInt) {
	  if (rndexp % 2) 
	       all [3] = new exprConst (-COUENNE_INFINITY);
	  else all [3] = new exprMin (new exprClone (all [1]),
				      new exprPow (new exprCopy (lbbase), 
						   new exprConst (expon)));
	}
	else all [3] = new exprCopy (all [1]);

	// all [5] is the lower bound value when lbbase <= ubbase <= 0

	if (expon > COUENNE_EPS) {
	  if (isInt && !(rndexp % 2))
	       all [5] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	  else all [5] = new exprConst (0);
	}
	else {
	  if (isInt || isInvInt) {
	    if (rndexp % 2)
	         all [5] = new exprPow (new exprClone (ubbase), new exprConst (expon));
	    else all [5] = new exprPow (new exprClone (lbbase), new exprConst (expon));
	  }
	  else all [5] = new exprConst (0);
	}

	lb = new exprMin (all, 6);

	// And now the upper bound ///////////////////////////////////

	if (expon > 0) {

	  // special case: bounds are referred to bounds only

	  ub = new exprMax (new exprPow (new exprClone (lbbase), new exprConst (expon)),
			    new exprPow (new exprClone (ubbase), new exprConst (expon)));
	} else {

	  expression **alu = new expression * [6];

	  alu [0] = new exprClone (all [0]);
	  alu [2] = new exprConst (0);
	  alu [4] = new exprClone (ubbase);

	  if ((expon > COUENNE_EPS) || ((isInt || isInvInt) && !(rndexp % 2)))
	    alu [1] = new exprPow (new exprCopy (ubbase), new exprConst (expon));
	  else alu [1] = new exprPow (new exprCopy (lbbase), new exprConst (expon));

	  // alu [3] is upper bound when lbbase <= 0 <= ubbase

	  if (expon < - COUENNE_EPS) 
	    alu [3] = new exprConst (COUENNE_INFINITY);
	  else if (isInt && !(rndexp % 2))
	    alu [3] = new exprPow (new exprMax (new exprCopy (lbbase), new exprCopy (ubbase)),
				   new exprConst (expon));
	  else alu [3] = new exprPow (new exprCopy (ubbase), new exprConst (expon));

	  // alu [5] is the upper bound value when lbbase <= ubbase <= 0

	  if (expon > COUENNE_EPS) {

	    if (isInt && !(rndexp % 2)) 
	      alu [5] = new exprPow (new exprClone(ubbase), new exprConst(expon));
	    else alu [5] = new exprConst (0);
	  }
	  else {
	    if (isInt || isInvInt) 
	      alu [5] = new exprPow (new exprClone(ubbase), new exprConst(expon));
	    else alu [5] = new exprConst (COUENNE_INFINITY);
	  }

	  ub = new exprMin (alu, 6);
	}
      }
    }
  }
}
