/*
 * Name:    exprLog.h
 * Author:  Pietro Belotti
 * Purpose: definition of logarithm
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRLOG_H
#define COUENNE_EXPRLOG_H

#include <exprInv.h>
#include <expression.h>
#include <math.h>


// class logarithm

class exprLog: public exprUnary {

 public:

  // Constructors, destructor
  exprLog  (expression *al): 
    exprUnary (al, log) {} //< non-leaf expression, with argument list

  ~exprLog () {}

  void print (std::ostream&);

  // differentiation
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // return expression of this same type with argument arg
  inline expression *mirror (expression *arg)
    {return new exprLog (arg);}

  // return derivative of univariate function of this type
  inline expression *mirror_d (expression *arg)
    {return new exprInv (arg);}

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
