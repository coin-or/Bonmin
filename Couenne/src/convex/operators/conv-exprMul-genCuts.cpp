/*
 * Name:    conv-exprMul-genCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: method to convexify multiplications
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprMul.h>
#include <exprDiv.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

/// add tangent at intersection with bounding box
void addImplTangent (const CouenneCutGenerator *, OsiCuts &, 
		     CouNumber, CouNumber, CouNumber, int, int, int, int);

/// generate convexification cut for constraint w = x*y

void exprMul::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // get bounds of numerator and denominator

  expression *xle, *xue, 
             *yle, *yue, 
             *wle, *wue;

  expression *xe = arglist_ [0];
  expression *ye = arglist_ [1];

  int wi = w  -> Index (), 
      xi = xe -> Index (), 
      yi = ye -> Index ();

  // if expression is x*c or c*y, with c constant from the problem
  // definition or from the branching rules, the expression is
  // linear. Add one convexification equality constraint.

  // check if either operand is constant

  bool is0const = (xe -> Type () == CONST),
       is1const = (ye -> Type () == CONST);

  CouNumber c0, c1;

  // is one of the two constant?

  if (is0const) c0 = xe -> Value ();  
  if (is1const) c1 = ye -> Value ();  

  // compute bounds

  xe -> getBounds (xle, xue);
  ye -> getBounds (yle, yue);
  w  -> getBounds (wle, wue);

  CouNumber xl = (*xle) (), xu = (*xue) (), 
            yl = (*yle) (), yu = (*yue) (),
            wl = (*wle) (), wu = (*wue) ();

  delete xle; delete xue;
  delete yle; delete yue;
  delete wle; delete wue;

  // check if either operator got constant because of the branching
  // rules: 

  bool i0s, i1s = i0s = false;

  if (!is0const) {

    if (is1const) i0s = (fabs (c1) * (xu-xl) < COUENNE_EPS);
    else          i0s = ((yu-yl)   * (xu-xl) < COUENNE_EPS);

    if (i0s) 
      c0 = 0.5 * (xl+xu);
  }

  // and y

  if (!is1const) {

    if (is0const) i1s = (fabs (c0) * (yu-yl) < COUENNE_EPS);
    else          i1s = ((xu-xl)   * (yu-yl) < COUENNE_EPS);

    if (i1s) 
      c1 = 0.5 * (yl+yu);
  }

  if (i0s) is0const = true;
  if (i1s) is1const = true;

  // right now c0 and c1 only have a value if the corresponding
  // expression is constant

  if (is0const || is1const) {

    if (cg -> isFirst () ||            // if first call or
	((xe -> Type () != CONST) &&   // neither term is a defined constant
	 (ye -> Type () != CONST))) {  // (=> implied by branching rule)

      if (is0const && is1const)

	// w = c0*c1, which is either because the intervals got very
	// narrow, or because these are indeed constant (which should
	// have been dealt with in simplify(), but who knows...)

	cg -> createCut (cs, c0 * c1, 0, wi, 1.);

      else {

	CouNumber coe;
	int ind;

	if (is0const) {coe = c0; ind = yi;} // c*y
	else          {coe = c1; ind = xi;} // x*c

	cg -> createCut (cs, 0., 0, wi, -1., ind, coe);
      }
    }

    return;
  }

  // Add McCormick convexification cuts:
  //
  // 1) w >= yl x + xl y - yl xl
  // 2) w >= yu x + xu y - yu xu
  //
  // 3) w <= yl x + xu y - yl xu
  // 4) w <= yu x + xl y - yu xl
  //
  // These cuts are added if the corresponding bounds are finite

  if (is_boundbox_regular (yl, xl)) cg -> createCut (cs, yl*xl, -1, wi, -1., xi, yl, yi, xl);
  if (is_boundbox_regular (yu, xu)) cg -> createCut (cs, yu*xu, -1, wi, -1., xi, yu, yi, xu);
  if (is_boundbox_regular (yl, xu)) cg -> createCut (cs, yl*xu, +1, wi, -1., xi, yl, yi, xu);
  if (is_boundbox_regular (yu, xl)) cg -> createCut (cs, yu*xl, +1, wi, -1., xi, yu, yi, xl);

  // extra cuts for the two cases w = xy >= l > 0 and w = xy <= u < 0
  //
  // These cuts are added at the tangent of the curve xy=k, on the
  // intersection with the bounding box. These are in fact NOT implied
  // by the above cuts (as happens for division, for instance) and may
  // be of help.

  // TODO: are these really useful/correct?

  if (wu < - COUENNE_EPS) {
    // check points A and B: second orthant intersections
    if ((xu*yl > wu) && (xl*yu <= wu)) {
      addImplTangent (cg, cs, xl, yl, wu, xi, yi, +1, +1); // A
      addImplTangent (cg, cs, xu, yu, wu, xi, yi, -1, +1); // B
    }
    else
      // check points C and D: fourth orthant intersections
      if ((xl*yu > wu) && (xu*yl <= wu)) {
	addImplTangent (cg, cs, xl, yl, wu, xi, yi, -1, -1); // C
	addImplTangent (cg, cs, xu, yu, wu, xi, yi, +1, -1); // D
      }
  }
  else 
    if (wl > COUENNE_EPS) {
      // check points A and B: third orthant intersections
      if ((xl*yl >= wl) && (xu*yu < wl)) {
	addImplTangent (cg, cs, xl, yu, wl, xi, yi, -1, -1); // A
	addImplTangent (cg, cs, xu, yl, wl, xi, yi, +1, -1); // B
      }
      else
	// check points C and D: first orthant intersections
	if ((xu*yu >= wl) && (xl*yl < wl)) {
	  addImplTangent (cg, cs, xl, yu, wl, xi, yi, +1, +1); // C
	  addImplTangent (cg, cs, xu, yl, wl, xi, yi, -1, +1); // D
	}
    }
}


/// add tangent to feasibility region of w=x*y at intersection with bounding box

void addImplTangent (const CouenneCutGenerator *cg, OsiCuts &cs,
		     CouNumber xb,    // x bound
		     CouNumber yb,    // y 
		     CouNumber wb,    // w 
		     int xi, int yi,  // indices of variables passed to createCut
		     int sign_check,  // get maximum or minimum of abscissa
		     int sign_ineq) { // lower or upper half-plane

  CouNumber xp, yp; // intersection point

  //  point is (xb, wb/xb) if xb*yb > wb, (wb/yb,yb) otherwise

  CouNumber check = xb*yb - wb; // violation (signed)

  if (sign_check < 0) // intersection point depends on violation check
    check = -check;

  if (check > 0) {xp = xb;    yp = wb/xb;}
  else           {xp = wb/yb; yp = yb;}

  // infinite or null bounds give vertical or horizontal cuts, useless
  if ((fabs (xp) < COUENNE_EPS) ||
      (fabs (yp) < COUENNE_EPS) ||
      (fabs (xp) > COUENNE_INFINITY) ||
      (fabs (yp) > COUENNE_INFINITY))
    return;

  CouNumber w_xp = wb / xp;

  cg -> createCut (cs, yp+w_xp, sign_ineq, yi, 1., xi, w_xp/xp);
}
