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
    exprUnary (al) {} //< non-leaf expression, with argument list
  ~exprOpp () {}

  // cloning method
  expression *clone () const
    {return new exprOpp (argument_ -> clone ());}

  // String equivalent (for comparisons)
  const std::string name () const {return exprUnary::name ("opp");}

  /// the operator's function
  inline unary_function F () {return opp;}

  // I/O
  void print (std::ostream&) const;

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

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
