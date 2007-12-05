/*
 * Name:    drawCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: print function and cut information to be later displayed
 *          through PlotUtils' graph
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "exprAux.hpp"


// [cool!] print graph-readable output for displaying inequalities
// on a Cartesian plane

void draw_cuts (OsiCuts &cs, const CouenneCutGenerator *cg, int j, expression *w, expression *img) {

  static bool first_draw = true;
  static CouNumber maxY = -COUENNE_INFINITY,
                   minY =  COUENNE_INFINITY;

  if ((img -> code () == COU_EXPRSIN) || 
      (img -> code () == COU_EXPRPOW) || 
      (img -> code () == COU_EXPREXP) || 
      (img -> code () == COU_EXPRLOG) || 
      (img -> code () == COU_EXPRCOS)) {

    //    fprintf (stderr, " ==> "); w -> print (std::cerr); fprintf (stderr, "\n");

    expression *lbe, *ube;

    expression *indep = img -> Argument ();

    if (!indep) 
      indep = img -> getFixVar ();

    int xi = indep -> Index ();

    //    fprintf (stderr, "looking into w_%d = f (x_%d)\n", 
    //	     w -> Index (), xi);

    indep -> getBounds (lbe, ube);

    CouNumber lb = (*lbe) (),
      ub = (*ube) ();

    delete lbe;
    delete ube;

    if (xi >= 0) {

      CouNumber curx = cg -> Problem () -> X (xi);

#define N_STEPS 100

      // plot function

      if (first_draw) {

	first_draw = false;

	for (CouNumber x = lb; 
	     x <= ub + COUENNE_EPS; 
	     x += ((ub - lb) / N_STEPS)) {

	  cg -> Problem () -> X () [xi] = x;
	  
	  CouNumber y = (*img) ();

	  if (y > maxY) maxY = y;
	  if (y < minY) minY = y;

	  fprintf (stderr, "%.12e %.12e\n", x, y);
	}

	maxY += (maxY-minY) / 20;
	minY -= (maxY-minY) / 20;
      }
	
      lb -= (ub-lb) / 20;
      ub += (ub-lb) / 20;

      // plot lines defining constraint (only for cuts involving at
      // most two variables (that is, w is a unary function)
      for (int jj=j; jj < cs.sizeRowCuts (); jj++) {

	CouNumber lb0 = lb, 
	  ub0 = ub;

	const double *el  = cs.rowCutPtr (jj) -> row (). getElements ();
	double  rhs = cs.rowCutPtr (jj) -> rhs ();

	if (fabs (el [1]) > COUENNE_EPS) {
	  lb0 = CoinMax (lb, CoinMin ((rhs - el[0] * minY) / el [1], (rhs - el[0] * maxY) / el [1]));
	  ub0 = CoinMin (ub, CoinMax ((rhs - el[0] * minY) / el [1], (rhs - el[0] * maxY) / el [1]));
	}

	fprintf (stderr, "#m=2,S=%d\n", (cs.rowCutPtr (jj) -> sense () == 'L') ? 10:11);

	fprintf (stderr, "%.12e %.12e\n", lb0, (rhs - el [1] * lb0) / el [0]);
	fprintf (stderr, "%.12e %.12e\n", ub0, (rhs - el [1] * ub0) / el [0]);
      }

      cg -> Problem () -> X () [xi] = curx;
      exit(0);
    }
  }
}
