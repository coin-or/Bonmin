/*
 * Name:    branchExprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for quadratic forms
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "exprQuad.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"


//#define DEBUG

/// set up branching object by evaluating many branching points for
/// each expression's arguments
CouNumber exprQuad::selectBranch (expression *w, 
				  const OsiBranchingInformation *info,
				  int &ind, 
				  double * &brpts, 
				  int &way) {

  const CouNumber *x = info -> solution_, 
                  *l = info -> lower_,
                  *u = info -> upper_; 

  CouNumber delta    = (*w) () - (*this) (), 
           *alpha    = (delta < 0) ? dCoeffLo_ : dCoeffUp_,
            maxcontr = -COUENNE_INFINITY;

  int bestVar = 
    nqterms_ ? *qindexI_ : 
    nlterms_ ? *index_   : -1; // initialize it to something very default

  if (!alpha) { // no convexification available,
                // find the largest interval

    CouNumber maxcontr = -COUENNE_INFINITY; // maximum contribution to current value

    int bestVar = -1;

    if (dIndex_) { // use indices in dIndex_ -- TODO: dIndex_ should
		   // be created by constructor

      for (int i=nDiag_; i--;) {

	int ind = dIndex_ [i];
	CouNumber diff;

#ifdef DEBUG
	printf ("[%10g %10g %10g] %4d %10g\n", x [ind], l [ind], u [ind], bestVar, maxcontr);
#endif

	//if ((diff = CoinMin (xi - l [qi], u [qi] - xi)) > maxcontr) {bestVar = qi; maxcontr = diff;}
	if ((diff = u [ind] - l [ind])                  > maxcontr) {bestVar = ind; maxcontr = diff;}
      }
    } else
      for (int i=nqterms_; i--;) {

	// find largest (in abs. value) coefficient with at least one
	// index within bounds

	int qi = qindexI_ [i],
            qj = qindexJ_ [i];

	CouNumber 
	  xi = x [qi], 
	  xj = x [qj], diff;

	if ((diff = CoinMin (xi - l [qi], u [qi] - xi)) > maxcontr) {bestVar = qi; maxcontr = diff;}
	if ((diff = CoinMin (xj - l [qj], u [qj] - xj)) > maxcontr) {bestVar = qj; maxcontr = diff;}
      }

    ind = bestVar;
    brpts = (double *) realloc (brpts, sizeof (double));
    //    *brpts = x [bestVar];
    *brpts = (l [bestVar] + u [bestVar]) / 2;
    way = TWO_RAND;

#ifdef DEBUG
    printf ("brExprQuad (%s): |delta| = %g, brpt = %g (%g), var = x%d, way = %d -- NO alpha\n",
	    (delta < 0) ? "lower" : "upper",
	    fabs (delta), *brpts, x [bestVar], bestVar, way);
#endif

    return fabs (delta);
  }

  int *indices = dIndex_;

  // there is a convexification already, find i = argmin {alpha_i (x_i
  // - l_i) (u_i - x_i)}

  for (int i=0; i < nDiag_; i++, indices++) {
    
    CouNumber curx = x [*indices],
      contrib = *alpha++ * 
      (curx         - l [*indices]) * 
      (u [*indices] - curx);

    if (fabs (contrib) > maxcontr) {
      bestVar  = *indices;
      maxcontr = contrib;
    }
  }

  ind = bestVar;
  way = TWO_RAND;

  brpts = (double *) realloc (brpts, sizeof (double));

  *brpts = midInterval (x [bestVar], l [bestVar], u [bestVar]);

  /*  if ((*brpts > ub - COUENNE_NEAR_BOUND) ||
      (*brpts < lb + COUENNE_NEAR_BOUND)) 

      *brpts = 0.5 * (lb + ub);*/

#ifdef DEBUG
  printf ("brExprQuad: |delta| = %g, brpt = %g (%g), var = x%d, way = %d\n",
	  fabs (delta), *brpts, x [bestVar], bestVar, way);
#endif

  return fabs (delta);
}
