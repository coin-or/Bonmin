/*
 * Name:    exprSub.h
 * Author:  Pietro Belotti
 * Purpose: definition of subtractions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUB_H
#define COUENNE_EXPRSUB_H

#include <exprOp.h>
#include <CouennePrecisions.h>
#include <CouenneProblem.h>

// class for subtraction

class exprSub: public exprOp {

 public:

  // Constructors, destructor
  exprSub  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprSub (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  ~exprSub () {}

  void print (std::ostream&);

  // function for the evaluation of the expression
  CouNumber operator () ();

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
  virtual inline int Linearity () {

    int lin1 = arglist_ [0] -> Linearity ();
    int lin2 = arglist_ [1] -> Linearity ();

    if (lin1 < lin2) return lin2;
    else             return lin1;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};


// compute subtraction

inline CouNumber exprSub::operator () () {

  exprOp:: operator () ();

  register CouNumber ret = - *sp--;
  return (currValue_ = ret + *sp--);
}

#endif
