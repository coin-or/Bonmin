/*
 * Name:    conv-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: utility methods to convexify multiplications
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"
#include "exprMul.hpp"
#include "exprBMul.hpp"
#include "exprConst.hpp"
#include "exprPow.hpp"
#include "exprClone.hpp"
#include "CouenneProblem.hpp"


/// check if two arguments point to the same variable

inline bool areSameVariables (expression *v1, expression *v2) {

  register int t1 = v1 -> Type (), t2;
  return (((t1 == VAR) || (t1 == AUX)) &&
	  (((t2 = v2 -> Type ()) == VAR) || (t2 == AUX)) && 
	  (v1 -> Index () == v2 -> Index ()));
}


/// Create standard formulation of this expression

exprAux *exprMul::standardize (CouenneProblem *p, bool addAux) {

  exprOp::standardize (p);

  if (nargs_ == 1)  // TODO: what happens really?
    return NULL;
  /* {
     exprAux *aux = arglist_ [0];
     arglist_ [0] = NULL;
     return aux;
     } */

  //expression *aux = arglist_ [0]; // why not this one?
  expression *aux = new exprClone (arglist_ [0]);

  for (int i = 1; i < nargs_ - 1; i++)
    aux = (areSameVariables (aux, arglist_ [i])) ? 
      (p -> addAuxiliary (new exprPow (new exprClone (aux), new exprConst (2.)))) : 
      (p -> addAuxiliary (new exprMul (new exprClone (aux), new exprClone (arglist_ [i]))));

  if (areSameVariables (aux, arglist_ [nargs_ - 1])) 
    aux    = new exprPow (new exprClone (aux), new exprConst (2.));
  else aux = new exprMul (new exprClone (aux), new exprClone (arglist_ [nargs_ - 1]));

  return (addAux ? (p -> addAuxiliary (aux)) : new exprAux (this, p -> domain ()));
}


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (expression *&lb, expression *&ub) {

  int i;

  if ((arglist_ [i=0] -> Type () == CONST) ||
      (arglist_ [i=1] -> Type () == CONST)) {

    CouNumber c = arglist_ [i] -> Value ();

    if (!i && (arglist_ [1] -> Type () == CONST)) { 

      // !i means i==0, or the first is constant. If you are here,
      // both are constant, which should not happen...

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


/// get lower/upper bounds of product f(x) g(x) in expression form

void exprMul::getBounds (CouNumber &lb, CouNumber &ub) {

  CouNumber lb1, ub1, lb2, ub2;

  arglist_ [0] -> getBounds (lb1, ub1);
  arglist_ [1] -> getBounds (lb2, ub2);

  if (ub1 < 0) { // use lb1, dominant
    if      (ub2 < 0) {lb = ub1*ub2; ub = lb1*lb2;}
    else if (lb2 > 0) {lb = lb1*ub2; ub = ub1*lb2;}
    else              {lb = lb1*ub2; ub = lb1*lb2;}
  } else if (lb1 > 0) { // use ub1, dominant
    if      (ub2 < 0) {lb = ub1*lb2; ub = lb1*ub2;}
    else if (lb2 > 0) {lb = lb1*lb2; ub = ub1*ub2;}
    else              {lb = ub1*lb2; ub = ub1*ub2;}
  } else { // there is a zero to consider
    if      (ub2 < 0) {lb = ub1*lb2; ub = lb1*lb2;}
    else if (lb2 > 0) {lb = lb1*ub2; ub = ub1*ub2;}
    else              {lb = CoinMin (lb1*ub2, lb2*ub1); 
                       ub = CoinMax (lb1*lb2, ub1*ub2);}
  }
}
