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

//#define NEW_TRIG

/// convex cuts for sine or cosine
void trigEnvelope (const CouenneCutGenerator *, OsiCuts &,
		   exprAux *, expression *, enum cou_trig);


// generate convexification cut for constraint w = sin (this)

void exprSin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

#ifdef NEW_TRIG
  trigEnvelope (cg, cs, w, w -> Image () -> Argument (), COU_SINE);
#else
  addHexagon (cg, cs, sin, w, w -> Image () -> Argument());
#endif
}


// generate convexification cut for constraint w = cos (this)

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

#ifdef NEW_TRIG
  trigEnvelope (cg, cs, w, w -> Image () -> Argument (), COU_COSINE);
#else
  addHexagon (cg, cs, cos, w, w -> Image () -> Argument());
#endif
}


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

  // add the lower envelope
  cg -> createCut (cs, f (lb) - lb, -1, w_ind, 1., x_ind, -1.); // left:  w - x <= f lb - lb 
  cg -> createCut (cs, f (ub) + ub, -1, w_ind, 1., x_ind,  1.); // right: w + x <= f ub + ub 

  // add the upper envelope
  cg -> createCut (cs, f (ub) - ub, +1, w_ind, 1., x_ind, -1.); // right: w - x >= cos ub - ub 
  cg -> createCut (cs, f (lb) + lb, +1, w_ind, 1., x_ind,  1.); // left:  w + x >= cos lb + lb 

  delete lbe;
  delete ube;
}


/// normalize angle within [0,b] (typically, pi or 2pi)
inline CouNumber modulo (register CouNumber a, register CouNumber b)
  {return a - b * floor (a/b);}


/// restrict to quarter of the interval [0,2pi]
int bayEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int, 
		 CouNumber, CouNumber, CouNumber, bool &, bool &);


/// real linearization of sine/cosine

void trigEnvelope (const CouenneCutGenerator *cg, // cut generator that has called us
		   OsiCuts &cs,                   // cut set to be enriched
		   exprAux *w,
		   expression *arg,
		   enum cou_trig which_trig) {

  expression *lbe, *ube;

  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) (),
            // if cosine, scale variables to pretend this is a sine problem
            displacement = (which_trig == COU_COSINE) ? M_PI_2 : 0.;

  delete lbe;
  delete ube;

  int xi = arg -> Index (),
      wi = w   -> Index ();

  // true if, in the first call (lb), a lower/upper chord was added
  // --> no such chord must be generated in the second call (ub)
  bool skip_up = false, 
       skip_dn = false;

  bayEnvelope (cg, cs, wi, xi, lb, ub, displacement, skip_up, skip_dn); // lb
  bayEnvelope (cg, cs, wi, xi, ub, lb, displacement, skip_up, skip_dn); // ub
}


// study single bay ( \__/ or /~~\ ) of the trigonometric function

int bayEnvelope (const CouenneCutGenerator *cg, // cut generator that has called us
		 OsiCuts &cs,                   // cut set to be enriched
		 int wi,
		 int xi,
		 CouNumber x0,           // starting point
		 CouNumber x1,           // other bound
		 CouNumber displacement, 
		 bool &skip_up, 
		 bool &skip_dn) {

  CouNumber tpt,
    rx0  = modulo (x0 + displacement, 2*M_PI),
    rx1  = rx0 + x1 - x0,
    base = x0 + displacement - rx0,
    sinrx0 = sin (rx0), zero;

  int
    up   = (modulo (rx0, 2*M_PI) < M_PI) ? +1 : -1,
    left = (x0 < x1)                     ? +1 : -1;

  zero = (up>0) ? 0. : M_PI;

  bool *s0, *s1;

  if (up>0) {s0 = &skip_up; s1 = &skip_dn;}
  else      {s0 = &skip_dn; s1 = &skip_up;}

  if (left * (modulo (rx0, M_PI) - M_PI_2) < 0) { 

    // after  flex (i.e., at \_ or /~ ) for left  bound, 
    // before flex (i.e., at _/ or ~\ ) for right bound

    // out of the "belly": tangent. If on upper bay we consider the
    // lower half-plane, and viceversa --> use -up
    cg -> addTangent (cs, wi, xi, x0, sin (rx0), cos (rx0), -up);

    // leftmost extreme to search for tangent point
    CouNumber extr0 = .75 * M_PI * (left+1) - M_PI_2 * up; 

    // in:
    if ((left * (rx1 - M_PI * ((left - up) / 2 + 1)) <= 0) ||   // if rx1 in same "belly", or
	(left * (rx1 - (tpt = trigNewton
			(rx0, extr0, extr0 + M_PI_2))) <= 0)) { // before closest leaning point 
      if (!*s1) // -> chord, if not already added in previous call
	*s1 = (cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1,         sin (rx1), up) > 0);
    } else     cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);
  }
  else {
    // after stationary point (i.e., _/ or ~\ ) for left bound, before
    // for right bound
  
    if (left * (rx1 - left * (zero + 5*M_PI_2)) < 0) {
      CouNumber cosrx0 = cos (rx0);
      if (up * (sin (rx1) - sinrx0 - cosrx0 * (rx1-rx0)) < 0) 
	// (b,sinb) below tangent --> tangent
	cg -> addTangent (cs, wi, xi, x0, sinrx0, cosrx0, up);
      else     // up: either chord or leaning plane
	if (left * (rx1 - (tpt = trigNewton (rx0, left * zero, left * (zero + M_PI_2)))) < 0) {
	  if (!*s0)
	    *s0 = cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1, sin (rx1), up) > 0;
	}
	else cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);
    } else   cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);

    // down: other chord or leaning plane
    if ((left * (rx1 - (zero + M_PI)) < 0) || 
	(left * (rx1 - (tpt = trigNewton (rx0, left * zero, left * (zero + M_PI_2)))) < 0)) {
      if (!*s1) 
	*s1 = (cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1,         sin (rx1), up) > 0);
    } else     cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);
  }

  return 0;
}
