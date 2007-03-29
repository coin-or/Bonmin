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


// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {
  return (((v1 -> Type () == VAR) || (v1 -> Type () == AUX)) &&
	  ((v2 -> Type () == VAR) || (v2 -> Type () == AUX)) && 
	  (v1 -> Index () == v2 -> Index ()));
}


// Create standard formulation of this expression

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


// get lower/upper bounds of product f(x) g(x) in expression form

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


// defined in conv-exprDiv.cpp
//bool is_boundbox_regular (CouNumber, CouNumber);


// generate convexification cut for constraint w = x*y

void exprMul::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {


  // TODO: add cuts considering w's lower and upper bounds (see
  // implied bounds)

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

  // if expression is x*c or c*y, with c constant from the problem
  // definition or from the branching rules, the expression is
  // linear. Add one convexification equality constraint.

  // check that either operand is constant

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

	if ((cut = cg -> createCut (c0 * c1, 0, w_ind, CouNumber (1.))))
	  cs.insert (cut);
      }
      else {

	CouNumber coe;
	int ind;

	if (is0const) {coe = c0; ind = y_ind;} // c*y
	else          {coe = c1; ind = x_ind;} // x*c

	if ((cut = cg -> createCut (CouNumber (0.), 0, w_ind, 
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
            yl = (*yle) (), yu = (*yue) ();

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
      && (cut = cg -> createCut (yl*xl, -1, w_ind, CouNumber (-1.), 
				 x_ind, yl, y_ind, xl)))
    cs.insert (cut);

  // 2)
  if (is_boundbox_regular (yu, xu)
      && (cut = cg -> createCut (yu*xu, -1, w_ind, CouNumber (-1.), 
				 x_ind, yu, y_ind, xu)))
    cs.insert (cut);

  // 3)
  if (is_boundbox_regular (yl, xu)
      && (cut = cg -> createCut (yl*xu, +1, w_ind, CouNumber (-1.), 
				 x_ind, yl, y_ind, xu)))
    cs.insert (cut);

  // 4)
  if (is_boundbox_regular (yu, xl)
      && (cut = cg -> createCut (yu*xl, +1, w_ind, CouNumber (-1.), 
				 x_ind, yu, y_ind, xl)))
    cs.insert (cut);

  delete xle; delete xue;
  delete yle; delete yue;
  delete wle; delete wue;
}
