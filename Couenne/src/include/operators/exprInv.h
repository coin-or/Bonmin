/*
 * Name:    exprInv.h
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRINV_H
#define COUENNE_EXPRINV_H

#include <exprUnary.h>


// the operator itself

inline CouNumber inv (register CouNumber arg) 
{return 1.0 / arg;}


// class inverse (1/f(x))

class exprInv: public exprUnary {

 public:

  // Constructors, destructor
  exprInv  (expression *al): 
    exprUnary (al, inv) {} //< non-leaf expression, with argument list
  ~exprInv () {}

  // cloning method
  expression *clone () const
    {return new exprInv (argument_ -> clone ());}

  // output "1/argument"
  void print (std::ostream&);

  // differentiation
  expression *differentiate (int index); 

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  virtual inline int Linearity () {
    if (argument_ -> Type () == CONST) return CONSTANT;
    else                               return NONLINEAR;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // construct linear over-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int upperLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
