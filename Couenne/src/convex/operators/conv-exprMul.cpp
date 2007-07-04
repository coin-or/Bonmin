/*
 * Name:    conv-exprMul.cpp
 * Author:  Pietro Belotti
 * Purpose: utility methods to convexify multiplications
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprMul.hpp>
#include <exprBMul.hpp>
#include <exprConst.hpp>
#include <exprPow.hpp>
#include <exprClone.hpp>
#include <CouenneProblem.hpp>


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

  int i;

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
