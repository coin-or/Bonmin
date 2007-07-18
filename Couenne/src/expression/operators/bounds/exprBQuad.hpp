/*
 * Name:    exprBQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of quadratic forms
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBQUAD_H
#define COUENNE_EXPRBQUAD_H

#include <exprOp.hpp>
#include <exprQuad.hpp>


/// method to actually compute the bound
CouNumber computeQBound (int, exprQuad *);

/// class to compute lower bound of a fraction based on the bounds of
/// both numerator and denominator

class exprLBQuad: public exprOp {

 public:

  /// Constructors, destructor
  exprLBQuad  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone () const
    {return new exprLBQuad (clonearglist (), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () () {return computeQBound (-1, NULL);}

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "LB_Quad";}
};


/// class to compute lower bound of a fraction based on the bounds of
/// both numerator and denominator

class exprUBQuad: public exprOp {

 public:

  /// Constructors, destructor
  exprUBQuad  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone () const
    {return new exprUBQuad (clonearglist (), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () () {return computeQBound (1, NULL);}

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "UB_Quad";}
};

#endif

