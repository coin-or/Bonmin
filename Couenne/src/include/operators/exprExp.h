/*
 * Name:    exprExp.h
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPREXP_H
#define COUENNE_EXPREXP_H

#include <exprUnary.h>
#include <math.h>

// class for the exponential

class exprExp: public exprUnary {

 public:

  // Constructors, destructor
  exprExp  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprExp (argument_ -> clone ());}

  /// the operator's function
  inline unary_function F () {return exp;}

  // output
  void print (std::ostream&) const;

  // differentiation
  expression *differentiate (int index); 

  // return expression of this same type with argument arg
  inline expression *mirror (expression *arg)
    {return new exprExp (arg);}

  // return derivative of univariate function of this type
  inline expression *mirror_d (expression *arg)
    {return new exprExp (arg);}

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  ///
  virtual enum expr_type code () {return COU_EXPREXP;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, char *);
};

#endif
