/*
 * Name:    conv-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to convexify multiplications
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprMul.h>
#include <exprBMul.h>
#include <exprConst.h>
#include <exprPow.h>
#include <exprDiv.h>
#include <exprClone.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


/// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {
  return (((v1 -> Type () == VAR) || (v1 -> Type () == AUX)) &&
	  ((v2 -> Type () == VAR) || (v2 -> Type () == AUX)) && 
	  (v1 -> Index () == v2 -> Index ()));
}


/// Create standard formulation of this expression

exprAux *exprMul::standardize (CouenneProblem *p) {

  exprOp::standardize (p);

  if (nargs_==1) return NULL;
  /* {
     exprAux *aux = arglist_ [0];
     arglist_ [0] = NULL;
     return aux;
     } */

  expression *aux = new exprClone (arglist_ [0]);

  for (int i=1; i < nargs_ - 1; i++)
    if (areSameVariables (aux, arglist_ [i]))
         aux = p -> addAuxiliary (new exprPow (aux, new exprConst (2)));
    else aux = p -> addAuxiliary (new exprMul (aux, new exprClone (arglist_ [i])));

  if (areSameVariables (aux, arglist_ [nargs_ - 1]))
       return  p -> addAuxiliary (new exprPow (aux, new exprConst (2)));
  else return  p -> addAuxiliary (new exprMul (aux, new exprClone (arglist_ [nargs_ - 1])));
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (expression *&lb, expression *&ub) {

  int i=2;

  if ((arglist_ [i=0] -> Type () == CONST) ||
      (arglist_ [i=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [i] -> Value ();

    if (!i && (arglist_ [1] -> Type () == CONST)) { 

      // !i means i==0, or the first is constant. If you are here,
      // both are constant, which should not happen. Anyway...

      CouNumber prod = c * arglist_ [1] -> Value ();

      lb = new exprConst (prod);
      ub = new exprConst (prod);

      return;
    }
    else {

      // expression is of the type c*x

      expression *lbi, *ubi;
      arglist_ [1-i] -> getBounds (lbi, ubi);

      if (c >= 0) {
	lb = new exprMul (new exprConst (c), lbi);
	ub = new exprMul (new exprConst (c), ubi);
      } else {
	lb = new exprMul (new exprConst (c), ubi);
	ub = new exprMul (new exprConst (c), lbi);
      }
    }
  }
  else {

    // expression is of the type x*y

    expression **almin = new expression * [4];
    expression **almax = new expression * [4];

    arglist_ [0] -> getBounds (almin [0], almin [1]);
    arglist_ [1] -> getBounds (almin [2], almin [3]);

    almax [0] = new exprClone (almin [0]);
    almax [1] = new exprClone (almin [1]);
    almax [2] = new exprClone (almin [2]);
    almax [3] = new exprClone (almin [3]);

    lb = new exprLBMul (almin, 4);
    ub = new exprUBMul (almax, 4);
  }
}


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

  OsiRowCut *cut;

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

  // check if either operator got constant because of the branching
  // rules: x

  if (is0const) c0 = xe -> Value ();  
  else {

    expression *xle, *xue;
    xe -> getBounds (xle, xue);

    c0 = (*xle) ();

    is0const = (fabs (c0 - (*xue) ()) < COUENNE_EPS);

    delete xle; delete xue;
  }

  // and y

  if (is1const) c1 = ye -> Value ();  
  else {

    expression *yle, *yue;
    ye -> getBounds (yle, yue);

    c1 = (*yle) ();

    is1const = (fabs (c1 - (*yue) ()) < COUENNE_EPS);

    delete yle; delete yue;
  }

  // right now c0 and c1 only have a value if the corresponding
  // expression is constant

  if (is0const || is1const) {

    if (cg -> isFirst () ||            // if first call or
	((xe -> Type () != CONST) &&   // neither term is a defined constant
	 (ye -> Type () != CONST))) {  // (and hence this follows from
				       // branching rule)

      if (is0const && is1const) {

	// strange case: w = c0*c1, should have been dealt with in
	// simplify, but who knows...

	if ((cut = cg -> createCut (c0 * c1, 0, wi, CouNumber (1.))))
	  cs.insert (cut);
      }
      else {

	CouNumber coe;
	int ind;

	if (is0const) {coe = c0; ind = yi;} // c*y
	else          {coe = c1; ind = xi;} // x*c

	if ((cut = cg -> createCut (CouNumber (0.), 0, wi, 
				    CouNumber (-1.), ind, coe)))
	  cs.insert (cut);
      }
    }

    return;
  }

  xe -> getBounds (xle, xue);
  ye -> getBounds (yle, yue);
  w  -> getBounds (wle, wue);

  CouNumber xl = (*xle) (), xu = (*xue) (), 
            yl = (*yle) (), yu = (*yue) (),
            wl = (*wle) (), wu = (*wue) ();

  delete xle; delete xue;
  delete yle; delete yue;
  delete wle; delete wue;

  // Add McCormick convexification cuts:
  //
  // 1) w >= yl x + xl y - yl xl
  // 2) w >= yu x + xu y - yu xu
  //
  // 3) w <= yl x + xu y - yl xu
  // 4) w <= yu x + xl y - yu xl
  //
  // if the corresponding bounds are finite

  // 1)
  if (is_boundbox_regular (yl, xl)
      && (cut = cg -> createCut (yl*xl, -1, wi, -1., xi, yl, yi, xl)))
    cs.insert (cut);

  // 2)
  if (is_boundbox_regular (yu, xu)
      && (cut = cg -> createCut (yu*xu, -1, wi, -1., xi, yu, yi, xu)))
    cs.insert (cut);

  // 3)
  if (is_boundbox_regular (yl, xu)
      && (cut = cg -> createCut (yl*xu, +1, wi, -1., xi, yl, yi, xu)))
    cs.insert (cut);

  // 4)
  if (is_boundbox_regular (yu, xl)
      && (cut = cg -> createCut (yu*xl, +1, wi, -1., xi, yu, yi, xl)))
    cs.insert (cut);

  // extra cuts for the two cases w = xy >= l > 0 and w = xy <= u < 0
  //
  // These cuts are added at the tangent of the curve xy=k, on the
  // intersection with the bounding box. These are in fact NOT implied
  // by the above cuts (as happens for division, for instance) and may
  // be of help.

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

  OsiRowCut *cut;
  CouNumber xp, yp; // intersection point

  //  point is (xb, wb/xb) if xb*yb > wb, (wb/yb,yb) otherwise

  CouNumber check = xb*yb - wb; // violation (signed)

  if (sign_check < 0) check = -check; // intersection point depends on
				      // the violation check

  if (check > 0) {xp = xb;    yp = wb/xb;}
  else           {xp = wb/yb; yp = yb;}

  // infinite or null bounds give vertical or horizontal cuts, useless
  if ((fabs (xp) < COUENNE_EPS) ||
      (fabs (yp) < COUENNE_EPS) ||
      (fabs (xp) > COUENNE_INFINITY-1) ||
      (fabs (yp) > COUENNE_INFINITY-1))
    return;

  CouNumber w_xp = wb / xp;

  if ((cut = cg -> createCut (yp+w_xp, sign_ineq, yi, 1., xi, w_xp/xp)))
    cs.insert (cut);
}
