/*
 * Name:    exprMul.h
 * Author:  Pietro Belotti
 * Purpose: definition of multiplications
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRMUL_H
#define COUENNE_EXPRMUL_H

#include <exprOp.h>
#include <exprAux.h>
#include <CouenneProblem.h>

// class for multiplications

class exprMul: public exprOp {

 public:

  // Constructors, destructor
  exprMul  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprMul (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  ~exprMul () {}

  // print expression
  void print (std::ostream&);

  // function for the evaluation of the expression
  inline CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is:
  virtual int Linearity ();

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  int lowerLinearHull (exprAux *, int *&, expression ***&, 
		       int **&, expression **&, enum con_sign *&);

  // construct linear over-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  int upperLinearHull (exprAux *, int *&, expression ***&, 
		       int **&, expression **&, enum con_sign *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};


// compute multiplication

inline CouNumber exprMul:: operator () () {

  exprOp:: operator () ();

  register CouNumber ret = *sp--;
  register int    n   =  nargs_;

  while (--n)
    ret *= *sp--;

  return (currValue_ = ret);
}

#endif
