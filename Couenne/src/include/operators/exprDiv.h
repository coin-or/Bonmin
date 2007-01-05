/*
 * Name:    exprDiv.h
 * Author:  Pietro Belotti
 * Purpose: definition of divisions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRDIV_H
#define COUENNE_EXPRDIV_H

#include <exprOp.h>
#include <CouennePrecisions.h>


// class for divisions

class exprDiv: public exprOp {

 public:

  // Constructors, destructor
  exprDiv (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprDiv (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  ~exprDiv () {}

  void print (std::ostream&);

  // function for the evaluation of the expression
  inline CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  inline int Linearity () {

    if (arglist_ [1] -> Type () == CONST)
         return arglist_ [0] -> Linearity ();
    else return NONLINEAR;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&lb, expression *&ub);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // construct linear over-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int upperLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};


// compute division

inline CouNumber exprDiv::operator () () {

  exprOp:: operator () ();

  register CouNumber denominator = *sp--;
  return (currValue_ = (*sp-- / denominator));
}

#endif
