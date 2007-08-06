/*
 * Name:    branchExprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch data for quadratic forms
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>
#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>


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

  CouNumber delta = (*w) () - (*this) (), 
    *alpha = (delta < 0) ? dCoeffLo_ : dCoeffUp_,
    maxcontr = -COUENNE_INFINITY;

  int bestVar = nqterms_ ? *qindexI_ : 
    nlterms_ ? *index_ : -1; // initialize it to something very default

  if (!alpha) { // no convexification available,
                // find the largest interval

    CouNumber maxcontr = -COUENNE_INFINITY; // maximum contribution to current value

    int bestVar = -1;

    for (int i=nqterms_; i--;) {

      // find largest (in abs. value) coefficient with at least one
      // index within bounds

      int qi = qindexI_ [i],
          qj = qindexJ_ [i];

      CouNumber xb = x [qi], diff;

      if ((diff = mymin (xb - l [qi], u [qi] - xb)) > maxcontr) {bestVar = qi; maxcontr = diff;}
      if ((diff = mymin (xb - l [qj], u [qj] - xb)) > maxcontr) {bestVar = qj; maxcontr = diff;}
    }

    ind = bestVar;
    brpts = (double *) realloc (brpts, sizeof (double));
    *brpts = x [bestVar];
    way = TWO_RAND;

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

  CouNumber lb = l [bestVar], 
            ub = u [bestVar];

  *brpts = x [bestVar];

  if ((*brpts > ub - COUENNE_NEAR_BOUND) ||
      (*brpts < lb + COUENNE_NEAR_BOUND)) 

    *brpts = 0.5 * (lb + ub);

  /*printf ("brExprQuad: |delta| = %g, brpt = %g (%g), var = x%d, way = %d\n",
    fabs (delta), *brpts, x [bestVar], bestVar, way);*/

  return fabs (delta);
}
