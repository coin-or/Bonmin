/*
 * Name:    conv-exprLog.C
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the logarithm operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprLog.h>
#include <exprConst.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// Lower convexification of exponential consists in (the north part
// of) the half-space defined by the segment from the point
// (lb,log(lb)) to the point (ub,log(ub))

int exprLog::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
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
  *sign           = COUENNE_GE;

  return 1;
}


// Upper convexification of logarithm is done by a set {1..nSamples()}
// of tangents to the log curve, computed by exprUnary::hull, each
// having three terms: auxiliary variable w, variable x which is the
// argument of exponential and whose coefficient is contained in
// coeffs [k] and whose rhs is contained in rhs [k], 0 <= k <
// nSamples().

int exprLog::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
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
    sign  [i] = COUENNE_LE;
  }

  delete [] x_coeff;

  return ns;
}


#define LOG_STEP 10
#define LOG_SCALE 1e-20

// generate convexification cut for constraint w = this

void exprLog::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
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

  int w_ind = aux       -> Index ();
  int x_ind = argument_ -> Index ();

  // fix lower bound

  if (l < COUENNE_EPS) 
    l = COUENNE_EPS;
  else   // lower segment (only put if lower bound is far enough from
	 // zero and upper is finite
    if ((u < COUENNE_INFINITY - 1) && 
	(check || ((w-log(l) * (u-l) < (x-l) * log(u)*log(l) - COUENNE_EPS)))) {

      cut   = new OsiRowCut;
      coeff = new CouNumber [2];
      index = new int       [2];

      CouNumber dx   = u-l;
      CouNumber logu = log (u);
      CouNumber dw   = logu - log (l);

      coeff [0] =  dx; index [0] = w_ind;
      coeff [1] = -dw; index [1] = x_ind;

      cut -> setLb (dx*logu - u*dw);
      cut -> setRow (2, index, coeff);

      printf ("Log lower: "); cut -> print ();

      cs.insert (cut);
    }

  // add tangent points: first choose sampling points

  int ns = cg -> nSamples ();

  if (u > COUENNE_INFINITY - 1)
    u = l + (LOG_STEP << ns);

  if (l < 2 * COUENNE_EPS)
    l = LOG_SCALE;

  if (x <= COUENNE_EPS)
    l = LOG_SCALE;

  if ((cg -> ConvType () == UNIFORM_GRID) || cg -> isFirst ()) {

    // choose sampling points. If unbounded, re-bound using a rule of
    // thumb where each point is taken every log 2 from the finite bound
 
    // now add tangent at each sampling point

    CouNumber sample = l, 
              step   = pow (u/l, 1/ns);

    for (int i = 0; i <= ns; i++) {

      addTangent (cs, w_ind, x_ind, sample, log (sample), 1./sample, -1);
      sample *= step;
    }
  } else if (cg -> ConvType () == CURRENT_ONLY)
    addTangent (cs, w_ind, x_ind, x, log (x), 1./x, -1);
  else {

    CouNumber sample = x;

    addTangent (cs, w_ind, x_ind, x, log (x), 1./x, -1);

    for (int i = 0; i <= ns/2; i++) {

      sample -= (x-l) / ns;
      addTangent (cs, w_ind, x_ind, sample, log (sample), 1./sample, -1);
    }

    sample = x;

    for (int i = 0; i <= ns/2; i++) {

      sample += (u-x) / ns;
      addTangent (cs, w_ind, x_ind, sample, log (sample), 1./sample, -1);
    }
  }
}
