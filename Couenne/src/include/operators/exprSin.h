/*
 * Name:    exprSin.h
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSIN_H
#define COUENNE_EXPRSIN_H

#include <math.h>

#include <exprUnary.h>
#include <exprConst.h>


// class for sin f(x)

class exprSin: public exprUnary {

 public:

  // Constructors, destructor
  exprSin  (expression *al): 
    exprUnary (al, sin) {} //< non-leaf expression, with argument list
  ~exprSin () {}

  void print (std::ostream&);

  // differentiation
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  virtual inline void getBounds (expression *&lb, expression *&ub)
    {lb = new exprConst (-1); ub = new exprConst (1);}

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
