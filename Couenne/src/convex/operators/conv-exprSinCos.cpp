/*
 * Name:    conv-exprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for sines and cosines
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiSolverInterface.hpp>
#include <CouenneTypes.h>
#include <CouenneCutGenerator.h>
#include <exprSin.h>
#include <exprCos.h>
#include <exprAux.h>


// generate convexification cut for constraint w = sin (this)

void exprSin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  trigGenCuts (w, cs, cg, sin);
}


// generate convexification cut for constraint w = cos (this)

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  trigGenCuts (w, cs, cg, cos);
}


// unified function, works with cos and sin

void trigGenCuts (exprAux *w, OsiCuts &cs, 
		  const CouenneCutGenerator *cg, unary_function f) {

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


// add lateral edges of the hexagon providing 

void addHexagon (const CouenneCutGenerator *cg, // cut generator that has called us
		 OsiCuts &cs,       // cut set to be enriched
		 unary_function f,  // sine or cosine
		 exprAux *aux,      // auxiliary variable
		 expression *arg) { // argument of cos/sin (should be a variable)

  OsiRowCut *cut;

  expression *lbe, *ube;
  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) ();

  // add the lower envelope, left: w - x <= f lb - lb 

  int x_ind = arg -> Index ();
  int w_ind = aux -> Index ();

  if ((cut = cg -> createCut (f (lb) - lb, -1, w_ind, CouNumber (1.),
			      x_ind, CouNumber (-1.))))
    cs.insert (cut);

  // and right: w + x <= f ub + ub 

  if ((cut = cg -> createCut (f (ub) + ub, -1, w_ind, CouNumber (1.),
			      x_ind, CouNumber (1.))))
    cs.insert (cut);

  // add the lower envelope, right: w - x >= cos ub - ub 

  if ((cut = cg -> createCut (f (ub) - ub, +1, w_ind, CouNumber (1.),
			      x_ind, CouNumber (-1.))))
    cs.insert (cut);

  // and left: w + x >= cos lb + lb 

  if ((cut = cg -> createCut (f (lb) + lb, +1, w_ind, CouNumber (1.),
			      x_ind, CouNumber (1.))))
    cs.insert (cut);
}
