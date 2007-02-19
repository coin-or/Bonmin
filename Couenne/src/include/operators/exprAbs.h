/*
 * Name:    exprAbs.h
 * Author:  Pietro Belotti
 * Purpose: definition of the absolute value of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRABS_H
#define COUENNE_EXPRABS_H

#include <math.h>

#include <exprUnary.h>
#include <exprConst.h>


// class for abs f(x)

class exprAbs: public exprUnary {

 public:

  // Constructors, destructor
  exprAbs  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list
  ~exprAbs () {}

  /// the operator's function
  inline unary_function F () {return fabs;}

  // cloning method
  expression *clone () const
    {return new exprAbs (argument_ -> clone ());}

  // String equivalent (for comparisons)
  const std::string name() const 
    {return "abs" + exprUnary::name();}

  // I/O
  void print (std::ostream&) const;

  // differentiation
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
