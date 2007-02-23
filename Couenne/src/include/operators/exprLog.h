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
    exprUnary (al) {} //< non-leaf expression, with argument list

  //  ~exprLog () {}

  // cloning method
  expression *clone () const
    {return new exprLog (argument_ -> clone ());}

  // String equivalent (for comparisons)
  const std::string name () const {return exprUnary::name ("log");}

  /// the operator's function
  inline unary_function F () {return log;}

  // I/O
  void print (std::ostream&) const;

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

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
