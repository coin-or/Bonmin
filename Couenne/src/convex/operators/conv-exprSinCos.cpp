/*
 * Name:    conv-exprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: convexification methods for sines and cosines
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
# define M_PI_2 1.57079632679489661923
#endif

#include <OsiSolverInterface.hpp>
#include <CouenneTypes.hpp>
#include <CouenneCutGenerator.hpp>
#include <exprSin.hpp>
#include <exprCos.hpp>
#include <exprAux.hpp>

#define NEW_TRIG

/// convex cuts for sine or cosine
int trigEnvelope (const CouenneCutGenerator *, OsiCuts &,
		   exprAux *, expression *, enum cou_trig);


/// add four cuts with slope 1 and -1
int addHexagon (const CouenneCutGenerator *, // cut generator that has called us
		 OsiCuts &,      // cut set to be enriched
		 enum cou_trig,  // sine or cosine
		 exprAux *,      // auxiliary variable
		 expression *);  // argument of cos/sin (should be a variable)


/// generate convexification cut for constraint w = sin (this)

void exprSin::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {

  int wi = w -> Index ();

  if (chg && !(cg -> isFirst ()) && 
      (chg [wi].lower == UNCHANGED) && 
      (chg [wi].upper == UNCHANGED))
    return;

#ifdef NEW_TRIG
  if (trigEnvelope (cg, cs, w, w -> Image () -> Argument (), COU_SINE) == 0)
#else
  if (addHexagon (cg, cs, COU_SINE, w, w -> Image () -> Argument()) == 0)
#endif
    {

    }
}


/// generate convexification cut for constraint w = cos (this)

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, int wind, 
			    CouNumber lbw, CouNumber ubw) {

  int wi = w -> Index ();

  if (chg && !(cg -> isFirst ()) && 
      (chg [wi].lower == UNCHANGED) && 
      (chg [wi].upper == UNCHANGED))
    return;

#ifdef NEW_TRIG
  if (trigEnvelope (cg, cs, w, w -> Image () -> Argument (), COU_COSINE) == 0) 
#else
  if (addHexagon (cg, cs, COU_COSINE, w, w -> Image () -> Argument()) == 0)
#endif
    {

    }
}


/// add lateral edges of the hexagon providing 

int addHexagon (const CouenneCutGenerator *cg, // cut generator that has called us
		 OsiCuts &cs,       // cut set to be enriched
		 enum cou_trig tt,  // sine or cosine
		 exprAux *aux,      // auxiliary variable
		 expression *arg) { // argument of cos/sin (should be a variable)

  // AW 2007-06-11: The following doesn't compile with MSVC++ because
  // sin and cos are ambiguous
  //unary_function fn = (tt == COU_SINE) ? sin : cos;

  // retrieve argument bounds
  expression *lbe, *ube;
  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) ();

  delete lbe;
  delete ube;

  int ncuts = 0,
    x_ind = arg -> Index (),
    w_ind = aux -> Index ();

  if (fabs (ub - lb) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (ub+lb), f, fp;

    if (tt == COU_SINE) {f = sin (x0); fp =  cos (x0);}
    else                {f = cos (x0); fp = -sin (x0);}

    return cg -> createCut (cs, f - fp*x0, 0, w_ind, 1., x_ind, -fp);
  }

  // add  /    \ envelope
  //      \    /

  // left
  if (lb > -COUENNE_INFINITY) { // if not unbounded
    if (tt == COU_SINE) {
      ncuts += cg -> createCut (cs, sin (lb) - lb, -1, w_ind, 1., x_ind, -1.); // up:  w-x <= f lb - lb 
      ncuts += cg -> createCut (cs, sin (lb) + lb, +1, w_ind, 1., x_ind,  1.); // dn:  w+x >= f lb + lb 
    }
    else {
      ncuts += cg -> createCut (cs, cos (lb) - lb, -1, w_ind, 1., x_ind, -1.); // up:  w-x <= f lb - lb 
      ncuts += cg -> createCut (cs, cos (lb) + lb, +1, w_ind, 1., x_ind,  1.); // dn:  w+x >= f lb + lb 
    }
  }

  // right
  if (ub <  COUENNE_INFINITY) { // if not unbounded
    if (tt == COU_SINE) {
      ncuts += cg -> createCut (cs, sin (ub) - ub, +1, w_ind, 1., x_ind, -1.); // dn: w - x >= f ub - ub 
      ncuts += cg -> createCut (cs, sin (ub) + ub, -1, w_ind, 1., x_ind,  1.); // up: w + x <= f ub + ub 
    }
    else {
      ncuts += cg -> createCut (cs, cos (ub) - ub, +1, w_ind, 1., x_ind, -1.); // dn: w - x >= f ub - ub 
      ncuts += cg -> createCut (cs, cos (ub) + ub, -1, w_ind, 1., x_ind,  1.); // up: w + x <= f ub + ub 
    }
  }

  return ncuts;
}


/// normalize angle within [0,b] (typically, pi or 2pi)
inline CouNumber modulo (register CouNumber a, register CouNumber b)
  {return a - b * floor (a/b);}


/// restrict to quarter of the interval [0,2pi]
int bayEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int, 
		 CouNumber, CouNumber, CouNumber, bool &, bool &);


/// real linearization of sine/cosine

int trigEnvelope (const CouenneCutGenerator *cg, // cut generator that has called us
		   OsiCuts &cs,                   // cut set to be enriched
		   exprAux *w,
		   expression *arg,
		   enum cou_trig which_trig) {

  expression *lbe, *ube;

  arg -> getBounds (lbe, ube);

  CouNumber lb = (*lbe) (), 
            ub = (*ube) (),
            // if cosine, scale variables to pretend this is a sine problem
            displ = (which_trig == COU_COSINE) ? M_PI_2 : 0.;

  delete lbe;
  delete ube;

  int ncuts = 0,
    xi = arg -> Index (),
    wi = w   -> Index ();

  if (fabs (ub - lb) < COUENNE_EPS) {

    CouNumber x0 = 0.5 * (ub+lb), f, fp;

    if (which_trig == COU_SINE) {f = sin (x0); fp =  cos (x0);}
    else                        {f = cos (x0); fp = -sin (x0);}

    return cg -> createCut (cs, f - fp*x0, 0, wi, 1., xi, -fp);
  }

  // true if, in the first call (lb), a lower/upper chord was added
  // --> no such chord must be generated in the second call (ub)
  bool skip_up = false, 
       skip_dn = false;

  if (lb > -COUENNE_INFINITY) ncuts += bayEnvelope (cg, cs, wi, xi, lb, ub, displ, skip_up, skip_dn);
  if (ub <  COUENNE_INFINITY) ncuts += bayEnvelope (cg, cs, wi, xi, ub, lb, displ, skip_up, skip_dn);

  return ncuts;
}


//                             __
// study single bay ( \__/ or /  \ ) of the trigonometric function
//

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
    base = x0 - rx0,
    sinrx0 = sin (rx0), zero;

  int ncuts = 0,
    up   = (rx0 < M_PI) ? +1 : -1,
    left = (x0  < x1)   ? +1 : -1;

  // starting point of the current bay
  zero = (up>0) ? 0. : M_PI;

  bool *s0, *s1;

  if (up>0) {s0 = &skip_up; s1 = &skip_dn;}
  else      {s0 = &skip_dn; s1 = &skip_up;}

  if (left * (modulo (rx0, M_PI) - M_PI_2) < 0) { 

    // after  flex (i.e., at \_ or /~ ) for left  bound, 
    // before flex (i.e., at _/ or ~\ ) for right bound

    // out of the "belly": tangent. If on upper bay we consider the
    // lower half-plane, and viceversa --> use -up
    ncuts += cg -> addTangent (cs, wi, xi, x0, sin (rx0), cos (rx0), -up);

    // leftmost extreme to search for tangent point
    CouNumber extr0 = .75 * M_PI * (left+1) - M_PI_2 * up; 

    // in:
    if ((left * (rx1 - M_PI * ((left - up) / 2 + 1)) <= 0) ||   // if rx1 in same "belly", or
	(left * (rx1 - (tpt = trigNewton
			(rx0, extr0, extr0 + M_PI_2))) <= 0)) { // before closest leaning point 
      if (!*s1) // -> chord, if not already added in previous call
	*s1 = ((ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1,  sin (rx1), up)) > 0);
    } else      ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);
  }
  else {

    // after  stationary point (i.e., _/ or ~\ ) for left  bound, 
    // before stationary point (i.e., /~ or \_ ) for right bound
  
    //    if (left * (rx1 - left * (zero + 5*M_PI_2)) < 0) {
    if (left * (rx1 - (4*left - up + 2) * M_PI_2) < 0) {
      CouNumber cosrx0 = cos (rx0);
      if (up * (sin (rx1) - sinrx0 - cosrx0 * (rx1-rx0)) < 0) 
	// (b,sinb) below tangent --> tangent
	ncuts += cg -> addTangent (cs, wi, xi, x0, sinrx0, cosrx0, -up);
      else {    // up: either chord or leaning plane
	CouNumber searchpt = M_PI_2 * (2 + 3*left - up);
	tpt = trigNewton (rx0, searchpt, searchpt + left * M_PI_2);
	if (left * (rx1 - tpt) < 0) {
	  if (!*s0)
	    *s0 = ((ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1, sin (rx1), -up)) > 0);
	} else      ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), -up);
      }
    } else {
      CouNumber searchpt = M_PI_2 * (2 + 3*left - up);
      tpt = trigNewton (rx0, searchpt, searchpt + left * M_PI_2);
      ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), -up);
    }

    // down: other chord or leaning plane
    if ((left * (rx1 - (zero + M_PI)) < 0) || 
	(left * (rx1 - (tpt = trigNewton (rx0, (2 +   left - up) * M_PI_2, 
		  			       (2 + 2*left - up) * M_PI_2))) < 0)) {
      if (!*s1) 
	*s1 = ((ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), x1, sin (rx1), up)) > 0);
    } else      ncuts += cg -> addSegment (cs, wi, xi, x0, sin (rx0), base + tpt, sin (tpt), up);
  }

  return ncuts;
}
