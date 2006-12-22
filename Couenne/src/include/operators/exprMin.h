/*
 * Name:    exprMin.h
 * Author:  Pietro Belotti
 * Purpose: definition of $\f(x_{\argmin_{i\in I} y_i})$ 
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRMIN_H
#define COUENNE_EXPRMIN_H

#include <exprOp.h>


//  class max

class exprMin: public exprOp {

 public:

  // Constructors, destructor
  exprMin  (expression **al, int n): 
    exprOp (al, n) {}
 
  exprMin  (expression *el0, expression *el1):
    exprOp (new expression * [4], 4) {
    arglist_ [0] = el0; arglist_ [1] = new exprCopy (el0);
    arglist_ [2] = el1; arglist_ [3] = new exprCopy (el1);
  }

  ~exprMin () {}

  // I/O
  void print (std::ostream &out)
    {exprOp:: print (out, (char *) "min", PRE);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // differentiation
  inline expression *differentiate (int) 
    {return NULL;} 

  // simplification
  inline expression *simplify () 
    {return NULL;}

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  virtual inline int Linearity () 
    {return NONLINEAR;}

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual inline exprAux *standardize (CouenneProblem *)
    {return NULL;}

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};


// compute maximum

// TODO: method that does not use the stack, computes elements only
// and returns evaluation of best val pointer


inline CouNumber exprMin::operator () () {

  exprOp:: operator () ();

  register CouNumber best_val = *sp--; 
  register CouNumber best_el  = *sp--; 

  int n = nargs_ / 2;

  while (--n) {

    register CouNumber val = *sp--;
    register CouNumber el  = *sp--;

    if (el < best_el) {
      best_el  = el;
      best_val = val;
    }
  }

  return (currValue_ = best_val);
}

#endif
