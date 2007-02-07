/*
 * Name:    exprSum.h
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUM_H
#define COUENNE_EXPRSUM_H

#include <exprOp.h>


//  class sum 

class exprSum: public exprOp {

 public:

  // Constructors, destructor
  exprSum  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprSum (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  ~exprSum () {}

  // cloning method
  expression *clone () const
    {return new exprSum (clonearglist (), nargs_);}

  // String equivalent (for comparisons)
  const std::string name() const {return "+" + exprOp::name();}

  // I/O
  virtual void print (std::ostream &);

  // function for the evaluation of the expression
  CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is:
  int Linearity ();

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


// compute sum

inline CouNumber exprSum::operator () () {

  exprOp:: operator () ();

  register CouNumber ret = *sp--; 
  register int       n   = nargs_;

  while (--n)
    ret += *sp--;

  return (currValue_ = ret);
}

#endif
