/*
 * Name:    exprInv.h
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRINV_H
#define COUENNE_EXPRINV_H

#include <exprUnary.h>


// the operator itself

inline CouNumber inv (register CouNumber arg) 
{return 1.0 / arg;}


// class inverse (1/f(x))

class exprInv: public exprUnary {

 public:

  // Constructors, destructor
  exprInv  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list
  ~exprInv () {}

  // cloning method
  expression *clone () const
    {return new exprInv (argument_ -> clone ());}

  /// the operator's function
  inline unary_function F () {return inv;}

  // String equivalent (for comparisons)
  const std::string name() const {return "inv" + exprUnary::name();}

  // output "1/argument"
  void print (std::ostream&) const;

  // differentiation
  expression *differentiate (int index); 

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  virtual inline int Linearity () {
    if (argument_ -> Type () == CONST) return CONSTANT;
    else                               return NONLINEAR;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
