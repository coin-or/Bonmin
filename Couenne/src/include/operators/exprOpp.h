/*
 * Name:    exprOpp.h
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPROPP_H
#define COUENNE_EXPROPP_H

#include <exprUnary.h>


// operator opp: returns the opposite of a number

inline CouNumber opp (register CouNumber arg) 
{return - arg;}


// class opposite 

class exprOpp: public exprUnary {

 public:

  // Constructors, destructor
  exprOpp  (expression *al): 
    exprUnary (al, opp) {} //< non-leaf expression, with argument list
  ~exprOpp () {}

  // cloning method
  expression *clone () const
    {return new exprOpp (argument_ -> clone ());}

  // String equivalent (for comparisons)
  const std::string name() const {return "opp" + exprUnary::name();}

  // I/O
  void print (std::ostream&);

  // differentiation
  expression *differentiate (int index); 

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  inline int Linearity ()
    {return argument_ -> Linearity ();}

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  inline int upperLinearHull (exprAux *, int *&, expression ***&, 
  //			      int **&, expression **&, enum con_sign *&)
  //    {return 0;}

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
