/*
 * Name:    conv-exprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for sines and cosines
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include <OsiSolverInterface.hpp>
#include <CouenneTypes.h>
#include <CouenneCutGenerator.h>
#include <exprSin.h>
#include <exprCos.h>
#include <exprAux.h>


// generate convexification cut for constraint w = sin (this)

void exprSin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  addHexagon (cg, cs, sin, w, w->Image()->Argument());
  //  trigGenCuts (w, cs, cg, sin);
}


// generate convexification cut for constraint w = cos (this)

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  addHexagon (cg, cs, cos, w, w->Image()->Argument());
  //  trigGenCuts (w, cs, cg, cos);
}


// unified function, works with cos and sin
/*
void trigGenCuts (exprAux *w, OsiCuts &cs, 
		  const CouenneCutGenerator *cg, unary_function f) {

  // Not needed as the same is done by getBounds and enforced by bound
  // propagation


  if (cg -> isFirst ()) { // Do not consider values of auxiliary variables 

    // No convexification is available yet, therefore add the six
    // convexification cuts for sine/cosine functions

    OsiRowCut *cut;

    int w_ind = w -> Index ();

    // A sine and cosine (loose) convex envelope is a hexagon, where
    // the two horizontal edges are the constraints -1 <= w <= 1,
    // while the other four edges have slope -1 or 1 (two of each).

    // Draw the two horizontal edges first: upper, w <= 1

    if ((cut = cg -> createCut (-1., +1, w_ind, CouNumber (1.), true)))
      cs.insert (cut);

    // and lower, w >= -1

    if ((cut = cg -> createCut (+1., -1, w_ind, CouNumber (1.), true)))
      cs.insert (cut);
  }

  // Now add the lower/upper envelope
  addHexagon (cg, cs, f, w, w->Image()->Argument());
}
*/

// add lateral edges of the hexagon providing 

void addHexagon (const CouenneCutGenerator *cg, // cut generator that has called us
		 OsiCuts &cs,       // cut set to be enriched
		 unary_function f,  // sine or cosine
		 exprAux *aux,      // auxiliary variable
		 expression *arg) { // argument of cos/sin (should be a variable)

  expression *lbe, *ube;
  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) ();

  int x_ind = arg -> Index ();
  int w_ind = aux -> Index ();

  // add the lower envelope, left: w - x <= f lb - lb 
  cg -> createCut (cs, f (lb) - lb, -1, w_ind, 1., x_ind, -1.);

  // and right: w + x <= f ub + ub 
  cg -> createCut (cs, f (ub) + ub, -1, w_ind, 1., x_ind,  1.);

  // add the lower envelope, right: w - x >= cos ub - ub 
  cg -> createCut (cs, f (ub) - ub, +1, w_ind, 1., x_ind, -1.);

  // and left: w + x >= cos lb + lb 
  cg -> createCut (cs, f (lb) + lb, +1, w_ind, 1., x_ind,  1.);

  delete lbe;
  delete ube;
}


///

CouNumber modulo2Pi (register CouNumber &a) {

  register CouNumber pi2 = 2 * M_PI;

  return a - pi2 * floor (a/pi2);
}

/// real linearization of sine/cosine

void trigEnvelope (const CouenneCutGenerator *cg, // cut generator that has called us
		   OsiCuts &cs,                   // cut set to be enriched
		   exprAux *w,
		   expression *arg,
		   enum cou_trig which_trig) {

  expression *lbe, *ube;
  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) ();

  int x_ind = arg -> Index ();
  int w_ind = w   -> Index ();

  // if cosine, scale variables to pretend this is a sine problem

  if (which_trig == COU_COSINE) {

    lb += M_PI_2;
    ub += M_PI_2;
  }

  CouNumber rlb   = modulo2Pi (lb),
            rub   = modulo2Pi (ub),
            delta = ub-lb;

  // check four cases of lb: 

  if        (rlb < M_PI_2)   { //    increasing, nonnegative


  } else if (rlb < M_PI)     { // nonincreasing, nonnegative
  } else if (rlb < 3*M_PI_2) { // nonincreasing,    negative
  } else                     { //    increasing,    negative
  }


  // do the same for ub
  delete lbe;
  delete ube;
}
