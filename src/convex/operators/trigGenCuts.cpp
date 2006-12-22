/*
 * Name:    trigGenCuts.C
 * Author:  Pietro Belotti
 * Purpose: common convexification for sine and cosine
 *
 * (C) Pietro Belotti, Carnegie Mellon University.
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>

#include <exprCos.h>
#include <exprAux.h>

#include <CouenneCutGenerator.h>


// unified function, works with cos and sin

void trigGenCuts (exprAux *w, OsiCuts &cs, 
		  const CouenneCutGenerator *cg, unary_function f) {

  if (cg -> isFirst ()) { // Do not consider values of auxiliary variables 

    // No convexification is available yet, therefore add the six
    // convexification cuts for sine/cosine functions

    OsiRowCut *cut;
    CouNumber *coeff;
    int       *index;

    // A sine and cosine (loose) convex envelope is a hexagon, where
    // the two horizontal edges are the constraints -1 <= w <= 1,
    // while the other four edges have slope -1 or 1 (two of each).

    // Draw the two horizontal edges first: upper, w <= 1

    cut   = new OsiRowCut;
    coeff = new CouNumber [1];
    index = new int [1];

    cut -> setUb (1.);
    coeff [0] = 1;
    index [0] = w -> Index ();

    cut -> setRow (1, index, coeff);

    cs.insert (cut);

    // and lower, w >= -1

    cut   = new OsiRowCut;
    coeff = new CouNumber [1];
    index = new int [1];

    cut -> setLb (-1.);
    coeff [0] = 1;
    index [0] = w -> Index ();

    cut -> setRow (1, index, coeff);

    cs.insert (cut);

    // Now add the lower/upper envelope

    addHexagon (cs, f, false, w, w -> Image () -> Argument ());
  }
  else { // yes, auxiliary variables' value can be considered

    // for now just add more hexagon constraints
    addHexagon (cs, f, true, w, w -> Image () -> Argument ());

    /*
    switch (cg -> ConvType ()) {
    case UNIFORM_GRID: 
    case CURRENT_ONLY: break;
    case AROUND_CURPOINT: break;
    }
    */
  }
}


// add lateral edges of the hexagon providing 

void addHexagon (OsiCuts &cs,       // cut set to be enriched
		 unary_function f,  // sine or cosine
		 bool check,        // should violation be checked before adding cut?
		 exprAux *aux,      // auxiliary variable
		 expression *arg) { // argument of cos/sin (should be a variable)

  OsiRowCut *cut;
  CouNumber *coeff;
  int       *index;

  expression *lbe, *ube;
  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) ();

  CouNumber w, x;

  if (check) {
    w = (*aux) ();
    x = (*arg) ();
  }

  printf ("Trig slopes: "); cut -> print ();

  // add the lower envelope, left: w - x <= f lb - lb 

  if (!check || (w - x > f (lb) - lb + COUENNE_EPS)) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [2];
    index = new int [2];

    coeff [0] =  1.; index [0] = aux -> Index ();
    coeff [1] = -1.; index [1] = arg -> Index ();

    cut -> setUb (f (lb) - lb);
    cut -> setRow (2, index, coeff);

    cs.insert (cut);
  }

  // and right: w + x <= f ub + ub 

  if (!check || (w + x > f (ub) + ub + COUENNE_EPS)) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [2];
    index = new int [2];

    coeff [0] = 1.; index [0] = aux -> Index ();
    coeff [1] = 1.; index [1] = arg -> Index ();

    cut -> setUb (f (ub) + ub);
    cut -> setRow (2, index, coeff);

    cs.insert (cut);
  }

  // add the lower envelope, right: w - x >= cos ub - ub 

  if (!check || (w - x < f (ub) - ub - COUENNE_EPS)) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [2];
    index = new int [2];

    coeff [0] =  1.; index [0] = aux -> Index ();
    coeff [1] = -1.; index [1] = arg -> Index ();

    cut -> setUb (f (ub) - ub);
    cut -> setRow (2, index, coeff);

    cs.insert (cut);
  }

  // and left: w + x >= cos lb + lb 

  if (!check || (w + x < f (lb) + lb - COUENNE_EPS)) {

    cut   = new OsiRowCut;
    coeff = new CouNumber [2];
    index = new int [2];

    coeff [0] =  1.; index [0] = aux -> Index ();
    coeff [1] = -1.; index [1] = arg -> Index ();

    cut -> setUb (f (lb) + lb);
    cut -> setRow (2, index, coeff);

    cs.insert (cut);
  }
}
