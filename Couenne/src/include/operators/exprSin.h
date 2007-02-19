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
    exprUnary (al) {} //< non-leaf expression, with argument list
  ~exprSin () {}

  // cloning method
  expression *clone () const
    {return new exprSin (argument_ -> clone ());}

  // String equivalent (for comparisons)
  const std::string name () const {return exprUnary::name ("sin");}

  /// the operator's function
  inline unary_function F () {return sin;}

  // I/O
  void print (std::ostream&) const;

  // differentiation
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  virtual inline void getBounds (expression *&lb, expression *&ub)
    {lb = new exprConst (-1); ub = new exprConst (1);}

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
