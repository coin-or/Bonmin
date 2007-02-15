/*
 * Name:    exprUnary.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class for univariate functions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRUNARY_H
#define COUENNE_EXPRUNARY_H

#include <iostream>

#include <expression.h>
#include <CouenneTypes.h>

//
// univariate operator-type expression: requires single argument. All
// unary functions are derived from this base class, which has a lot
// of common methods that need not be re-implemented by any univariate
// class.
//

class exprUnary: public expression {

 protected:

  expression     *argument_; //< single argument taken by this expression
  unary_function  f_;        //< single variable function (e.g. sin, log...)

 public:

  /// node type
  virtual inline enum nodeType Type () 
    {return UNARY;}

  /// Constructor
  exprUnary  (expression *argument, unary_function f): 
    argument_ (argument),        //< non-leaf expression, with argument list
    f_        (f)         {}

  /// Destructor
  ~exprUnary () 
    {if (argument_) delete argument_;}

  /// return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
    {return argument_;}

  // string equivalent
  virtual const std::string name () const;

  // I/O
  virtual void print (std::ostream &, const std::string &, enum pos) const;

  // compute value of unary operator
  virtual inline CouNumber operator () ()
    {return (currValue_ = f_ ((*argument_) ()));}

  // dependence on variable set
  bool inline dependsOn (int *list, int n) 
    {return argument_ -> dependsOn (list, n);}

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is:
  //
  // 0: a constant
  // 1: linear
  // 2: quadratic
  // 3: nonlinear non-quadratic
  //
  // for general univariate functions, return nonlinear.
  virtual inline int Linearity ()
    {return NONLINEAR;}

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar () {return argument_;}
};

#endif
