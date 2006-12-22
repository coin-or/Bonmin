/*
 * Name:    conv-exprExp.C
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the exponential operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprAux.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

#define MULT_FACTOR 10


// Lower convexification of exponential is done by a set
// {1..nSamples()} of tangents to the exponential curve, computed by
// exprUnary::hull, each having three terms: auxiliary variable w,
// variable x which is the argument of exponential, whose coefficient
// is contained in coeffs [k] and whose rhs is contained in rhs [k], 0
// <= k < nSamples().

int exprExp::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  int ns = nSamples ();

  expression **x_coeff = new expression * [ns];

  nterms = new int [ns];

  for (int i=0; i<ns; i++) 
    nterms [i] = 2;

  allocateCon (ns, nterms, coeff, indices, rhs, sign);

  hull (x_coeff, rhs);

  for (int i=0; i<ns; i++) {
    coeff [i] [0] = new exprConst (1); indices [i] [0] = w         -> Index ();
    coeff [i] [1] = x_coeff [i];       indices [i] [1] = argument_ -> Index ();
    sign  [i]     = COUENNE_GE;
  }

  delete [] x_coeff;

  return ns;
}


// Upper convexification of exponential consists in (the south part
// of) the half-space defined by the segment from the point
// (lb,exp(lb)) to the point (ub,exp(ub))

int exprExp::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  coeff    = new expression ** [1];
  *coeff   = new expression * [2];
  indices  = new int * [1];
  *indices = new int [2];
  rhs      = new expression * [1];
  sign     = new enum con_sign [1];
  nterms   = new int [1];

  segment (coeff [0] [1], *rhs);

  *nterms = 2;
  coeff   [0] [0] = new exprConst (1);
  indices [0] [0] = w         -> Index ();
  indices [0] [1] = argument_ -> Index ();
  *sign           = COUENNE_LE;

  return 1;
}


// generate convexification cut for constraint w = this

void exprExp::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber w = (*aux) (),
            x = (*argument_) (),
            l = (*le) (),
            u = (*ue) ();

  bool check = cg -> isFirst () || !(cg -> addViolated ());

  CouNumber *coeff;
  int       *index;
  OsiRowCut *cut;

  // upper segment

  if ((u < COUENNE_INFINITY - 1) 
      && (check || ((w-exp(l) * (u-l) > (x-l) * exp(u)*exp(l) + COUENNE_EPS)))) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [2];
    index = new int       [2];

    CouNumber dx   = u-l;
    CouNumber expu = exp (u);
    CouNumber dw   = expu - exp (l);

    coeff [0] =  dx; index [0] = aux       -> Index ();
    coeff [1] = -dw; index [1] = argument_ -> Index ();

    cut -> setUb (dx*expu - u*dw);
    cut -> setRow (2, index, coeff);

    printf ("Exp upper: "); cut -> print ();

    cs.insert (cut);
  }

  // lower convexification: start with trivial envelope w >= 0

  if (check || (w < - COUENNE_EPS)) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [1];
    index = new int       [1];

    coeff [0] = 1; 
    index [0] = aux -> Index ();

    cut -> setLb (0);
    cut -> setRow (1, index, coeff);

    printf ("Exp trivial: "); cut -> print ();

    cs.insert (cut);
  }

  // add tangent points: first choose sampling points

  int ns = cg -> nSamples ();

  if ((cg -> ConvType () == UNIFORM_GRID) || cg -> isFirst ()) {

    // choose sampling points. If unbounded, re-bound using a rule of
    // thumb where each point is taken every log 2 from the finite bound
 
    if (l < - COUENNE_INFINITY + 1) {
      if (u > COUENNE_INFINITY - 1) {

	l = - log (MULT_FACTOR) * ns/2;
	u =   log (MULT_FACTOR) * ns/2;
      } else 
	l = u - log (2) * ns;

    } else if (u > COUENNE_INFINITY - 1) 
	u = l + log (2) * ns;

    // now add tangent at each sampling point

    CouNumber sample = l, step = (u-l) / ns;

    for (int i = 0; i <= ns; i++) {

      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, exp (sample), exp (sample), +1);
      sample += step;
    }
  } else if (cg -> ConvType ()== CURRENT_ONLY)
    addTangent (cs, aux -> Index (), argument_ -> Index (), 
		x, exp (x), exp (x), +1);
  else {

    CouNumber sample = x;

    addTangent (cs, aux -> Index (), argument_ -> Index (), 
		x, exp (x), exp (x), +1);

    for (int i = 0; i <= ns/2; i++) {

      sample -= (x-l) / ns;
      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, exp (sample), exp (sample), +1);
    }

    sample = x;

    for (int i = 0; i <= ns/2; i++) {

      sample += (u-x) / ns;
      addTangent (cs, aux -> Index (), argument_ -> Index (), 
		  sample, exp (sample), exp (sample), +1);
    }
  }
}
