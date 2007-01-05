/*
 * Name:    conv-exprUnary.C
 * Author:  Pietro Belotti
 * Purpose: methods to convexify n-ary operators
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprUnary.h>
#include <exprSum.h>
#include <exprSub.h>
#include <exprClone.h>
#include <exprOpp.h>
#include <exprDiv.h>
#include <exprMul.h>

#include <CouenneProblem.h>


// Create standard formulation of this expression, by:
//
// - creating auxiliary w variables and corresponding expressions
// - returning linear counterpart as new constraint (to replace 
//   current one)

exprAux *exprUnary::standardize (CouenneProblem *p) {

  exprAux *subst;

  if ((subst = argument_ -> standardize (p)))
    argument_ = subst;

  return p -> addAuxiliary (this);
}


// Method to create a linear, convex (concave) hull of a given
// function that is convex (concave)
/*
void exprUnary::hull (expression **& coeff, expression **& rhs) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  // now generate a set of underestimating constraints 

  int ns = nSamples ();

  coeff = new expression * [ns];
  rhs   = new expression * [ns];

  ns--;

  for (int i=0; i<=ns; i++) {

    CouNumber alpha = (CouNumber) i / ns;

    expression *sample;

    if           (!i)  sample = uba;   // x_k = ub
    else if (i == ns)  sample = lba;   // x_k = lb
    else               sample =        // x_k = alpha lb + (1-alpha) ub, with 0 < alpha < 1

			 // introduce new operator?
      new exprSum
      (new exprMul    (new exprConst   (alpha), new exprClone (lba)),
       new exprMul    (new exprConst (1-alpha), new exprClone (uba)));

    expression *f_sample  = mirror   (sample);                 // value of f at sample: f  (x_k)
    expression *fp_sample = mirror_d (new exprClone (sample)); // derivative at x_k:    f' (x_k)

    coeff [i] = new exprOpp (new exprClone (fp_sample));

    rhs [i] = new exprSub (f_sample, new exprMul (fp_sample, new exprClone (sample)));
  }
}


// Method to create a segment approximating from above (below) a
// function that is convex (concave)

void exprUnary::segment (expression*& coeff, expression*& rhs) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  // now generate one overestimating constraint

  // this is the slope of the segment over the curve exp (x) from lba
  // to uba: (exp(uba) - exp (lba)) / (uba - lba)
 
  expression *slope = new exprDiv (new exprSub (mirror        (lba), mirror        (uba)),
				   new exprSub (new exprClone (uba), new exprClone (lba)));
  coeff = slope;

  rhs = new exprDiv (new exprSub (new exprMul (mirror (lba), new exprClone (uba)),
				  new exprMul (mirror (uba), new exprClone (lba))),
		     new exprSub (new exprClone (uba), new exprClone (lba)));
}
*/

