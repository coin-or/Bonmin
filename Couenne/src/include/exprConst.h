/*
 * Name:    exprConst.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprConst
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCONST_H
#define COUENNE_EXPRCONST_H

#include <iostream>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprCopy.h>


// constant-type operator

class exprConst: public expression {

 public:

  // node type
  inline enum nodeType Type () 
    {return CONST;}

  // value of expression
  inline CouNumber Value () const 
    {return currValue_;}

  // Constructors, destructor
  exprConst  (CouNumber value) 
    {currValue_ = value;}
  ~exprConst () {}

  // I/O
  inline void print (std::ostream &out) 
    {out << currValue_;}

  // return constant's value
  inline CouNumber operator() () 
    {return currValue_;}

  // differentiation
  inline expression *differentiate (int) 
    {return new exprConst (0);}

  // dependence on variable set
  inline bool dependsOn (int *, int) 
    {return false;}

  // simplify
  inline expression *simplify () 
    {return NULL;}

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  inline int Linearity ()
    {return CONSTANT;}

  // Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) {
    lower = new exprCopy (this);
    upper = new exprCopy (this);
  }

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  int lowerLinearHull (exprAux *, int *&, expression ***&, 
		       int **&, expression **&, enum con_sign *&);

  // Create standard formulation of this expression
  inline exprAux *standardize (CouenneProblem *)
    {return NULL;}

  // generate convexification cut for constraint w = this
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
