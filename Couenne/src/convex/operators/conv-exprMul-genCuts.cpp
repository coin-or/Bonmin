/*
 * Name:    conv-exprMul-genCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: method to convexify multiplications
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprMul.h>
#include <exprPow.h>
#include <exprDiv.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

/// add tangent at intersection with bounding box
void addImplTangent (const CouenneCutGenerator *, OsiCuts &, 
		     CouNumber, CouNumber, CouNumber, int, int, int, int);

/// used below to specify k in curve y=k/x
static CouNumber multInv;

inline CouNumber kinv   (register CouNumber x)
{return multInv / x;}

inline CouNumber kinvp  (register CouNumber x)
{return -multInv / (x*x);}

inline CouNumber kinvpp (register CouNumber x)
{return 2*multInv / (x*x*x);}


/// Add cut around curve x*y=k 

void contourCut (const CouenneCutGenerator *cg,
		 OsiCuts &cs, 
		 CouNumber xp, CouNumber yp, // current point
		 CouNumber wb,               // bound on w
		 int sign,                   // is wb lower or upper?
		 CouNumber x0, CouNumber y0, // (allegedly) outside point
		 CouNumber x1, CouNumber y1, //             inside
		 int xi, int yi, int wi) {   // indices of the variables

  // Upper right corner of the bounding box of (x,y) is feasible,
  // the opposite corner is not, hence there is a cut violated by
  // (x0,y0).

  // If (x0,y0) is not in the same orthant as the contour in
  // question, move it in so that we can apply a Newton step to
  // find closest point on contour.

  int xsign = (x1 >= 0) ? 1 : -1, // define orthant where the "inside
      ysign = (y1 >= 0) ? 1 : -1; // point" lies

  if      (((xsign > 0) ? xp : -xp) <= COUENNE_EPS)
    if    (((ysign > 0) ? yp : -yp) <= COUENNE_EPS) { 

      // opposite orthant, put in the right one where constraint is violated
      xp = yp = sqrt (fabs (wb))/2; 
      if (xsign<0) xp = -xp;
      if (ysign<0) yp = -yp;
    }                                                // otherwise, must cross one axis only:
    else                                            {xp = sqrt (fabs(wb/yp)); if (xsign<0) xp=-xp;}//y
  else if (((ysign > 0) ? yp : -yp) <= COUENNE_EPS) {yp = sqrt (fabs(wb/xp)); if (ysign<0) yp=-yp;}//x

  multInv = wb;

  CouNumber 
    // tangent point closest to current point
    xt    = powNewton (xp, yp, kinv, kinvp, kinvpp),
    // coefficient of w in the lifted cut
    alpha = ((fabs (x1) < COUENNE_INFINITY) && 
	     (fabs (y1) < COUENNE_INFINITY)) ? 
       ((2*wb/xt - y1 - wb*x1 / (xt*xt)) / (x1*y1 - wb)) : 0;

  //  printf ("+++++ %d %d %d. [%c] xp (%g,%g) wb %g out(%g,%g) in(%g,%g) --> [%g,%g] alpha %g\n",
  //	         xi, yi, wi, (sign<0) ? '-' : '+', xp, yp, wb, x0, y0, x1, y1, xt, wb/xt, alpha);

  if (alpha != 0)
    cg     -> createCut (cs, alpha*wb + 2*wb/xt, sign, wi, alpha, yi, 1., xi, wb/(xt*xt));
  else  cg -> createCut (cs,            2*wb/xt, sign,            yi, 1., xi, wb/(xt*xt));
}


/// generate convexification cut for constraint w = x*y

void exprMul::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg) {

  // TODO: unify with exprDiv::generateCuts

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

  bool cLX,  cRX,  cLY,  cRY,  cLW,  cRW = 
       cLX = cRX = cLY = cRY = cLW = true;

  if (!(cg -> isFirst ()) && chg) {
    cLX = chg [xi].lower != UNCHANGED;  cRX = chg [xi].upper != UNCHANGED;
    cLY = chg [yi].lower != UNCHANGED;  cRY = chg [yi].upper != UNCHANGED;
    cLW = chg [wi].lower != UNCHANGED;  cRW = chg [wi].upper != UNCHANGED;
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

  if ((cLX || cLY) && is_boundbox_regular (yl, xl)) cg -> createCut (cs, yl*xl,-1,wi,-1.,xi,yl,yi,xl);
  if ((cRX || cRY) && is_boundbox_regular (yu, xu)) cg -> createCut (cs, yu*xu,-1,wi,-1.,xi,yu,yi,xu);
  if ((cRX || cLY) && is_boundbox_regular (yl, xu)) cg -> createCut (cs, yl*xu,+1,wi,-1.,xi,yl,yi,xu);
  if ((cLX || cRY) && is_boundbox_regular (yu, xl)) cg -> createCut (cs, yu*xl,+1,wi,-1.,xi,yu,yi,xl);

  // add different cuts, to cut out current point in bounding box but
  // out of the hyperbola's belly

  CouNumber x0 = (*(arglist_ [0])) (),
            y0 = (*(arglist_ [1])) ();

  // If w=xy and w >= l > 0 (resp. w <= u < 0) are "tight" bounds
  // (i.e. they are tighter than those obtained through propagation of
  // x and y's bounds), McCormick's convexification is not tight as
  // the surface has a curve contour at w=l (resp. w=u).
  //
  // McCormick rules induce a tangent to this contour at the bounds of
  // both variables, but it may be useful to add further cuts along
  // the contour to eliminate infeasible point (x0,y0,w0), which may
  // be in the convexification but out of the contour (on its "convex"
  // side, or "out of the belly").
  //
  // Suppose P (xt,l/xt) (resp. (xt,u/xt) is the point on the contour
  // closest to (x0,y0), found through a Newton method. The cut is
  // tangent to the contour in P and has the form
  //
  //        y - l/xt >= -l/(xt^2) (x-xt)   if xl*yl < l and xu*yu > l
  //        y - l/xt <= -l/(xt^2) (x-xt)   if xl*yl > l and xu*yu < l
  //
  // (resp. y - u/xt <= -u/(xt^2) (x-xt)   if xl*yu > u and xu*yl < u
  //        y - u/xt >= -u/(xt^2) (x-xt)   if xl*yu < u and xu*yl > u)
  //
  // These can be lifted to satisfy, at equality, the point
  // (xu,yu,wu=xu*yu) (resp. (xl,yl,wl=xl*yl)), where xl and xu are
  // lower and upper bound of x, etc.
  //
  //        alpha (w - l) + y - l/xt >= -l/(xt^2) (x-xt) ...
  //
  // where alpha is such that the relation holds at equality at the
  // point (xu,yu,xu*yu):
  //
  //    alpha = [-yu + l/xt - l/(xt^2)(xu-xt)] / (xu*yu - l)

  if ((x0 > xl + COUENNE_EPS) && (y0 > yl + COUENNE_EPS) &&
      (x0 < xu + COUENNE_EPS) && (y0 < yu + COUENNE_EPS)) {

    if (cLW && (wl > 0) && (x0*y0 < wl)) { // that is, if (x0,y0) is out of the contour

      CouNumber xyl = xl * yl;

      // first and third orthant
      if      ((xyl <  wl) && (xu*yu >=wl)) contourCut (cg,cs, x0,y0, wl, +1, xl,yl, xu,yu, xi,yi,wi);
      else if ((xyl >= wl) && (xu*yu < wl)) contourCut (cg,cs, x0,y0, wl, -1, xu,yu, xl,yl, xi,yi,wi);
    }

  // Similarly for w <= u < 0 

    if (cRW && (wu < 0) && (x0*y0 > wu)) { // that is, if (x0,y0) is out of the contour

      CouNumber xuyl = xl * yu;

      // second and fourth orthant
      if      ((xuyl > wu) && (xl*yu <=wu)) contourCut (cg,cs, x0,y0, wu, +1, xu,yl, xl,yu, xi,yi,wi);
      else if ((xuyl <=wu) && (xl*yu > wu)) contourCut (cg,cs, x0,y0, wu, -1, xl,yu, xu,yl, xi,yi,wi);
    }
  }
}
