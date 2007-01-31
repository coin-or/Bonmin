/*
 * Name:    conv-exprMul.C
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
#include <exprClone.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {
  return ((((v1 -> Type () == VAR) || (v1 -> Type () == AUX)) &&
	   ((v2 -> Type () == VAR) || (v2 -> Type () == AUX))) 
	  && (v1 -> Index () == v2 -> Index ()));
}


// Create standard formulation of this expression

exprAux *exprMul::standardize (CouenneProblem *p) {

  exprOp::standardize (p);

  expression *aux = new exprClone (arglist_ [0]);

  for (int i=1; i < nargs_ - 1; i++)
    if (areSameVariables (aux, arglist_ [i]))
         aux = p -> addAuxiliary (new exprPow (aux, new exprConst (2)));
    else aux = p -> addAuxiliary (new exprMul (aux, new exprClone (arglist_ [i])));

  if (areSameVariables (aux, arglist_ [nargs_ - 1]))
       return  p -> addAuxiliary (new exprPow (aux, new exprConst (2)));
  else return  p -> addAuxiliary (new exprMul (aux, new exprClone (arglist_ [nargs_ - 1])));
}


// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (expression *&lb, expression *&ub) {

  int i;

  if ((arglist_ [i=0] -> Type () == CONST) ||
      (arglist_ [i=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [i] -> Value ();

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
  else {

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


// construct linear under-estimator for expression within problem *p
// (p is used to add convexification constraints)
/*
int exprMul::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // expression w = xy is convexified from below with two linear
  // constraints (lx and ly are lower bounds of x and y, ux and uy are
  // upper bounds
  //
  // w >= ly x + lx y - ly lx
  // w >= uy x + ux y - uy ux

  if ((arglist_ [0] -> Type () == CONST) ||
      (arglist_ [1] -> Type () == CONST))
    return 0;

  expression *lx, *ux, *ly, *uy;

  arglist_ [0] -> getBounds (lx, ux);
  arglist_ [1] -> getBounds (ly, uy);

  nterms = new int [2];
  nterms [0] = nterms [1] = 3;

  allocateCon (2, nterms, coeff, indices, rhs, sign);

  coeff [0] [0] = new exprConst (-1); indices [0] [0] = w            -> Index ();
  coeff [0] [1] = new exprClone (ly); indices [0] [1] = arglist_ [0] -> Index ();
  coeff [0] [2] = new exprClone (lx); indices [0] [2] = arglist_ [1] -> Index ();
  rhs   [0] = new exprMul (lx, ly);
  sign  [0] = COUENNE_LE;

  coeff [1] [0] = new exprConst (-1); indices [1] [0] = w            -> Index ();
  coeff [1] [1] = new exprClone (uy); indices [1] [1] = arglist_ [0] -> Index ();
  coeff [1] [2] = new exprClone (ux); indices [1] [2] = arglist_ [1] -> Index ();
  rhs   [1] = new exprMul (ux, uy);
  sign  [1] = COUENNE_LE;

  return 2;
}


// similarly, construct linear over-estimator for expression within
// problem *p (p is used to add convexification constraints). It is
// also used when this function appears with a minus sign in the
// expression

int exprMul::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // expression w = xy is convexified from below with two linear
  // constraints (lx and ly are lower bounds of x and y, ux and uy are
  // upper bounds
  //
  // w <= ly x + ux y - ly ux
  // w <= uy x + lx y - uy lx

  if ((arglist_ [0] -> Type () == CONST) ||
      (arglist_ [1] -> Type () == CONST))
    return 0;

  expression *lx, *ux, *ly, *uy;

  arglist_ [0] -> getBounds (lx, ux);
  arglist_ [1] -> getBounds (ly, uy);

  nterms = new int [2];
  nterms [0] = nterms [1] = 3;

  allocateCon (2, nterms, coeff, indices, rhs, sign);

  coeff [0] [0] = new exprConst (-1); indices [0] [0] = w            -> Index ();
  coeff [0] [1] = new exprClone (ly); indices [0] [1] = arglist_ [0] -> Index ();
  coeff [0] [2] = new exprClone (ux); indices [0] [2] = arglist_ [1] -> Index ();
  rhs   [0] = new exprMul (ly, ux);
  sign  [0] = COUENNE_GE;

  coeff [1] [0] = new exprConst (-1); indices [1] [0] = w            -> Index ();
  coeff [1] [1] = new exprClone (uy); indices [1] [1] = arglist_ [0] -> Index ();
  coeff [1] [2] = new exprClone (lx); indices [1] [2] = arglist_ [1] -> Index ();
  rhs   [1] = new exprMul (uy, lx);
  sign  [1] = COUENNE_GE;

  return 2;
}
*/

// is x finite?
bool is_finite (CouNumber x);

// generate convexification cut for constraint w = this

void exprMul::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  // get bounds of numerator and denominator

  expression *xle, *xue, 
             *yle, *yue, 
             *wle, *wue;

  expression *xe = arglist_ [0];
  expression *ye = arglist_ [1];

  OsiRowCut *cut;

  int w_ind = w  -> Index (), 
      x_ind = xe -> Index (), 
      y_ind = ye -> Index ();

  // if expression is c*x, with c constant, everything gets easier...

  if ((xe -> Type () == CONST) || 
      (ye -> Type () == CONST)) {

    if (cg -> isFirst ()) {

      if ((xe -> Type () == CONST) && 
	  (ye -> Type () == CONST)) {

	// strange case: w = c1*c2, should have been dealt with in
	// simplify, but who knows...

	if ((cut = cg -> createCut (xe -> Value () * ye -> Value (), 
				    0, w_ind, CouNumber (1.))))
	  cs.insert (cut);
      }
      else {

	CouNumber coe;
	int ind;

	if (xe -> Type () != CONST) {coe = ye -> Value (); ind = x_ind;}
	else                        {coe = xe -> Value (); ind = y_ind;}
	/*
	  printf ("========> w=cx: c=%.4f, x_%d, w_%d... ", coe, ind, w_ind);
	  xe -> print (std::cout); printf (" * ");
	  ye -> print (std::cout); printf ("\n");
	*/
	if ((cut = cg -> createCut (CouNumber (0.), 0, w_ind, CouNumber (-1.), ind, coe))){
	  //	  cut -> print ();
	  cs.insert (cut);
	}
      }
    }

    return;
  }

  xe -> getBounds (xle, xue);
  ye -> getBounds (yle, yue);
  w  -> getBounds (wle, wue);

  CouNumber xl = (*xle) (), xu = (*xue) (), 
            yl = (*yle) (), yu = (*yue) ();

  /*
  printf ("Mult: x = ");
  xe  -> print (std::cout); printf (" [");
  xle -> print (std::cout); printf (",");
  xue -> print (std::cout); printf ("]    y = ");

  ye  -> print (std::cout); printf (" [");
  yle -> print (std::cout); printf (",");
  yue -> print (std::cout); printf ("];   w = ");

  w   -> print (std::cout); printf ("\n");
  */

  // Add McCormick convexification cuts:
  //
  // 1) w >= yl x + xl y - yl xl
  // 2) w >= yu x + xu y - yu xu
  //
  // 3) w <= yl x + xu y - yl xu
  // 4) w <= yu x + xl y - yu xl


  // 1) 

  if (is_finite (yl) && is_finite (xl) 
      && (cut = cg -> createCut (yl*xl, -1, w_ind, CouNumber (-1.), x_ind, yl, y_ind, xl))) {

    //    printf ("--- cut 1: "); cut -> print ();
    cs.insert (cut);
  }

  // 2) 

  if (is_finite (yu) && is_finite (xu) 
      && (cut = cg -> createCut (yu*xu, -1, w_ind, CouNumber (-1.), x_ind, yu, y_ind, xu))) {

    //    printf ("--- cut 2: "); cut -> print ();
    cs.insert (cut);
  }

  // 3) 

  if (is_finite (yl) && is_finite (xu) 
      && (cut = cg -> createCut (yl*xu, +1, w_ind, CouNumber (-1.), x_ind, yl, y_ind, xu))) {

    //    printf ("--- cut 3: "); cut -> print ();
    cs.insert (cut);
  }

  // 4) 

  if (is_finite (yu) && is_finite (xl) 
      && (cut = cg -> createCut (yu*xl, +1, w_ind, CouNumber (-1.), x_ind, yu, y_ind, xl))) {

    //    printf ("--- cut 4: "); cut -> print ();
    cs.insert (cut);
  }
}
