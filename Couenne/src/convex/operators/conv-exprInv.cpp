/*
 * Name:    conv-exprInv.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the inverse operator
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"

#include "exprInv.hpp"
#include "exprClone.hpp"
#include "exprConst.hpp"
#include "exprMin.hpp"
#include "exprOpp.hpp"
#include "exprDiv.hpp"
#include "exprSum.hpp"
#include "exprMul.hpp"
#include "exprPow.hpp"

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

// compute upper- and lower bounds of the expression w = 1/f(x) given
// bounds of f(x)

void exprInv::getBounds (expression *&lb, expression *&ub) {

  expression *lba, *uba;
  argument_ -> getBounds (lba, uba);

  expression **all = new expression * [6];
  all [0] = new exprConst (0.);      all [1] = new exprConst (- COUENNE_INFINITY);  // l<0<u
  all [2] = new exprOpp   (lba);     all [3] = new exprInv   (uba);                 // 0<l<u
  all [4] = new exprClone (uba);     all [5] = new exprInv   (new exprClone (uba)); // l<u<0

  lb = new exprMin (all, 6);

  expression **alu = new expression * [6];
  alu [0] = new exprConst (0.);      alu [1] = new exprConst (COUENNE_INFINITY);   // l<0<u
  alu [2] = new exprClone (all [2]); alu [3] = new exprInv (new exprClone (lba));  // 0<l<u
  alu [4] = new exprClone (uba);     alu [5] = new exprInv (new exprClone (lba));  // l<u<0

  ub = new exprMin (alu, 6);
}


// compute VALUE of lower and upper bound of expression
void exprInv::getBounds (CouNumber &lb, CouNumber &ub) {

  register CouNumber lba, uba;

  argument_ -> getBounds (lba, uba);

  if ((uba < 0) || (lba > 0)) {
    lb = 1./uba;
    ub = 1./lba;
  } else {
    lb = -COUENNE_INFINITY;
    ub =  COUENNE_INFINITY;
  }
}


#define MIN_DENOMINATOR 1e-6

// generate convexification cut for constraint w = 1/x

void exprInv::generateCuts (expression *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {
  CouNumber l, u;
  argument_ -> getBounds (l, u);

  if ((l < - COUENNE_EPS) && (u > COUENNE_EPS)) // there is no convexification
    return;

  int wi = aux       -> Index (), 
      xi = argument_ -> Index ();

  bool cL = !chg || (cg -> isFirst ()) || (chg [xi].lower() != t_chg_bounds::UNCHANGED);
  bool cR = !chg || (cg -> isFirst ()) || (chg [xi].upper() != t_chg_bounds::UNCHANGED);

  // special case: l and u are very close, replace function with
  // linear term

  if (fabs (u - l) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (u+l);
    if (cL || cR) 
      cg -> createCut (cs, 2/x0, 0, wi, 1., xi, 1/(x0*x0));
    return;
  }

  // upper segment (or lower if x<0)

  if (cL || cR) {
    // bounding box is within ]0,+inf[
    if ((l> COUENNE_EPS) && (u< COU_MAX_COEFF)) cg -> createCut (cs, 1/l+1/u, -1, wi,1., xi,1/(l*u));
    if ((u<-COUENNE_EPS) && (l>-COU_MAX_COEFF)) cg -> createCut (cs, 1/l+1/u, +1, wi,1., xi,1/(l*u));
    // bounding box is within ]-inf,0[
  }

  // choose sampling points. 

  // if unbounded, use a rule of thumb
  int ns = cg -> nSamples ();
  if      (l < - COUENNE_INFINITY) l = ns * (u-1); // (-infinity, u] where u < 0
  else if (u >   COUENNE_INFINITY) u = ns * (l+1); // [l, +infinity) where l > 0

  // make bounds nonzero
  if (fabs (l) < COUENNE_EPS) l = (l<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;
  if (fabs (u) < COUENNE_EPS) u = (u<0) ? - MIN_DENOMINATOR : MIN_DENOMINATOR;

  // bound
  cg -> addEnvelope 
    (cs, (l > 0) ? +1 : -1, 
     inv, oppInvSqr, wi, xi, 
     (cg -> isFirst ()) ? // is this first call?
       // place it somewhere in the interval (we don't care)
       ((l > COUENNE_EPS) ? l : u) :
       // otherwise, replace it where it gives deepest cut
       powNewton ((*argument_) (), (*aux) (), inv, oppInvSqr, inv_dblprime),
     l, u, chg);
}
