/*
 * Name:    conv-exprDiv.C
 * Author:  Pietro Belotti
 * Purpose: standardization and convexification methods for divisions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOp.h>
#include <exprDiv.h>
#include <exprClone.h>
#include <exprMul.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

// Create standard formulation of this expression
exprAux *exprDiv::standardize (CouenneProblem *p) {

  exprOp::standardize (p);
  return p -> addAuxiliary (this);
}


// The equation q = n/d is given to be convexified. After computing
// the lower/upper bounds for q (which are infinite if the bounds of d
// are different in sign), we translate the equation into n = qd and
// apply the convexification for products.

int exprDiv::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // expression w = x/y is convexified by transforming it into wy = x
  // and applying the same convexification as for products

  expression *q, *Q, *d, *D;

  getBounds (q, Q); //

  arglist_ [1] -> getBounds (d, D);

  nterms = new int [2];
  nterms [0] = nterms [1] = 3;

  allocateCon (2, nterms, coeff, indices, rhs, sign);

  coeff [0] [0] = new exprConst (-1); indices [0] [0] = arglist_ [0] -> Index ();
  coeff [0] [1] = new exprClone (d);  indices [0] [1] = w            -> Index ();
  coeff [0] [2] = new exprClone (q);  indices [0] [2] = arglist_ [1] -> Index ();
  rhs   [0] = new exprMul (d, q);
  sign  [0] = COUENNE_LE;

  coeff [1] [0] = new exprConst (-1); indices [1] [0] = arglist_ [0] -> Index ();
  coeff [1] [1] = new exprClone (D);  indices [1] [1] = w            -> Index ();
  coeff [1] [2] = new exprClone (Q);  indices [1] [2] = arglist_ [1] -> Index ();
  rhs   [1] = new exprMul (D, Q);
  sign  [1] = COUENNE_LE;

  return 2;
}


// construct linear under-estimator for expression within problem *p
// (p is used to add convexification constraints)

int exprDiv::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // expression w = x/y is convexified by transforming it into wy = x
  // and applying the same convexification as for products

  expression *q, *Q, *d, *D;

  getBounds (q, Q); //

  arglist_ [1] -> getBounds (d, D);

  nterms = new int [2];
  nterms [0] = nterms [1] = 3;

  allocateCon (2, nterms, coeff, indices, rhs, sign);

  coeff [0] [0] = new exprConst (-1); indices [0] [0] = arglist_ [0] -> Index ();
  coeff [0] [1] = new exprClone (d);  indices [0] [1] = w            -> Index ();
  coeff [0] [2] = new exprClone (Q);  indices [0] [2] = arglist_ [1] -> Index ();
  rhs  [0] = new exprMul (d, Q);
  sign [0] = COUENNE_GE;

  coeff [1] [0] = new exprConst (-1); indices [1] [0] = arglist_ [0] -> Index ();
  coeff [1] [1] = new exprClone (D);  indices [1] [1] = w            -> Index ();
  coeff [1] [2] = new exprClone (q);  indices [1] [2] = arglist_ [1] -> Index ();
  rhs  [1] = new exprMul (D, q);
  sign [1] = COUENNE_GE;

  return 2;
}


// generate convexification cut for constraint w = x/y

void exprDiv::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // get bounds of numerator and denominator

  expression *yle, *yue;

  arglist_ [1] -> getBounds (yle, yue);

  CouNumber yl = (*yle) (), 
            yu = (*yue) ();

  // if the denominator's bound interval has 0 as internal point,
  // there is no convexification

  if ((yl < - COUENNE_EPS) && (yu > COUENNE_EPS)) 
    return;

  expression *xle, *xue, *wle, *wue;

  arglist_ [0] -> getBounds (xle, xue);
  w            -> getBounds (wle, wue);

  expression *xe = arglist_ [0];
  expression *ye = arglist_ [1];

  CouNumber x = (*xe) (), xl = (*xle) (), xu = (*xue) (),
            y = (*ye) (),
            w0 = (*w) (), wl = (*wle) (), wu = (*wue) ();

  // Add McCormick convexification cuts. Reduce w = x/y to the case
  // x = wy and apply the same rule as for multiplications:
  //
  // 1) x >= yl w + wl y - yl wl
  // 2) x >= yu w + wu y - yu wu
  //
  // 3) x <= yl w + wu y - yl wu
  // 4) x <= yu w + wl y - yu wl

  CouNumber *coeff;
  int       *index;
  OsiRowCut *cut;

  bool check = cg -> isFirst () || !(cg -> addViolated ());

  // 1) 

  if (check || (x < yl*w0 + wl*y - yl*wl - COUENNE_EPS)) {

    coeff = new CouNumber [3];
    index = new int       [3];
    cut   = new OsiRowCut;

    coeff [0] = -1; index [0] = xe -> Index ();
    coeff [1] = yl; index [1] = w  -> Index ();
    coeff [2] = wl; index [2] = ye -> Index ();

    cut -> setUb (yl * wl);
    cut -> setRow (3, index, coeff);

    printf ("Division1: "); cut -> print ();
  
    cs.insert (cut);
  }

  // 2) 

  if (check || (x < yu*w0 + wu*y - yu*wu - COUENNE_EPS)) {

    coeff = new CouNumber [3];
    index = new int       [3];
    cut   = new OsiRowCut;

    coeff [0] = -1; index [0] = xe -> Index ();
    coeff [1] = yu; index [1] = w  -> Index ();
    coeff [2] = wu; index [2] = ye -> Index ();

    cut -> setUb (yu * wu);
    cut -> setRow (3, index, coeff);

    printf ("Division2: "); cut -> print ();
    cs.insert (cut);
  }

  // 3) 

  if (check || (x > yl*w0 + wu*y - yl*wu + COUENNE_EPS)) {

    coeff = new CouNumber [3];
    index = new int       [3];
    cut   = new OsiRowCut;

    coeff [0] = -1; index [0] = xe -> Index ();
    coeff [1] = yl; index [1] = w  -> Index ();
    coeff [2] = wu; index [2] = ye -> Index ();

    cut -> setLb (yl * wu);
    cut -> setRow (3, index, coeff);

    printf ("Division3: "); cut -> print ();

    cs.insert (cut);
  }

  // 4) 

  if (check || (x < yu*w0 + wl*y - yu*wl - COUENNE_EPS)) {

    coeff = new CouNumber [3];
    index = new int       [3];
    cut   = new OsiRowCut;

    coeff [0] = -1; index [0] = xe -> Index ();
    coeff [1] = yu; index [1] = w  -> Index ();
    coeff [2] = wl; index [2] = ye -> Index ();

    cut -> setLb (yu * wl);
    cut -> setRow (3, index, coeff);

    printf ("Division4: "); cut -> print ();

    cs.insert (cut);
  }
}
