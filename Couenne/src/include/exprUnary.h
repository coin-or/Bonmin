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

  // node type
  virtual inline enum nodeType Type () 
    {return UNARY;}

  // Constructors, destructor
  exprUnary  (expression *argument, unary_function f): 
    argument_ (argument),        //< non-leaf expression, with argument list
    f_        (f)         {}

  ~exprUnary () 
    {if (argument_) delete argument_;}
  /*
  // return a pointer to an object with the same type of operator f_
  // as this class, and with arg as argument_
  virtual inline expression *mirror (expression *arg)
    {return new exprUnary (arg, f_);}

  // return a pointer to an object with the same type as the
  // derivative of operator f_ as this class, and with arg as argument_
  virtual inline expression *mirror_d (expression *arg)
    {return new exprUnary (arg, f_);}
  */
  // return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
    {return argument_;}

  // I/O
  virtual void print (std::ostream &, const std::string &, enum pos);

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

  // for auxiliary variable w = f(x), returns coefficient of x and
  // right hand side of the expression 
  //
  // w - (f(ub)-f(lb))/(ub-lb) * x  >=< f(lb) - (f(ub)-f(lb))/(ub-lb) * lb
  //
  // used as a convexification (lower bound) of a concave function
  // (and vice versa, as a concavification -- upper bound -- of a
  // convex function). The sign depends on the convexity of the
  // operator.
  void segment (expression *&, expression *&);

  // for auxiliary variable w = f(x), returns coefficients of x and
  // right hand sides of the expressions
  //
  // w - f'(x_k) * x  >=< f(x_k) - f'(x_k) * x_k
  //
  // where k is the set of sample points used to create convex
  // (concave) hull. Used as a convexification (lower bound) of a
  // concave function (and vice versa, as a concavification -- upper
  // bound -- of a convex function). The sign depends on the convexity
  // of the operator.
  void hull (expression **&, expression **&);
};

#endif
