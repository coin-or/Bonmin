/*
 * Name:    exprBound.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class for variable bounds 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBOUND_HPP
#define COUENNE_EXPRBOUND_HPP

#include <iostream>

#include <CouenneTypes.hpp>
#include <exprVar.hpp>

/// These are bound expression classes. They are used in the parametric
/// convexification part to obtain lower/upper bounds of an expression
/// as a function of the expression itself.
///
/// For example, the lower and upper bounds to expression (x1 - exp
/// (x2)) are (l1 - exp (u2)) and (u1 - exp (l2)), respectively, where
/// l1 (l2) is the lower bound of x1 (x2) and u1 (u2) is the upper
/// bound of x1 (x2).
///
/// A lower/upper bound of an expression is a function of all bounds in
/// the expression and is known only when all variables bounds are
/// known.


/// lower bound 

class exprLowerBound: public exprVar {

 public:

  /// Node type
  inline enum nodeType Type () 
    {return CONST;}

  /// Constructor
  exprLowerBound (int varIndex): 
    exprVar (varIndex) {}

  /// Copy constructor
  exprLowerBound (const exprLowerBound &src): 
    exprVar (src) {}

  /// cloning method
  exprLowerBound *clone () const
    {return new exprLowerBound (*this);}

  /// Print to iostream
  void print (std::ostream &out = std::cout, 
	      bool = false, CouenneProblem * = NULL) const
    {out << "l_" << varIndex_;}

  /// return the value of the variable
  inline CouNumber operator () () 
    {return (currValue_ = expression::lbounds_ [varIndex_]);}

  /// differentiation
  inline expression *differentiate (int) 
    {return new exprConst (0);}

  /// dependence on variable set
  inline int dependsOn (int *, int) 
    {return 0;}

  /// get a measure of "how linear" the expression is:
  virtual inline int Linearity () 
    {return CONST;}

  /// code for comparisons
  virtual enum expr_type code ()
    {return COU_EXPRLBOUND;}
};


/// upper bound 

class exprUpperBound: public exprVar {

 public:

  /// Node type
  inline enum nodeType Type () 
    {return CONST;}

  /// Constructor
  exprUpperBound (int varIndex): 
    exprVar (varIndex) {}

  /// Copy constructor
  exprUpperBound (const exprUpperBound &src): 
    exprVar (src) {}

  /// cloning method
  exprUpperBound *clone () const
    {return new exprUpperBound (*this);}

  /// Print to iostream
  void print (std::ostream &out = std::cout, 
	      bool = false, CouenneProblem * = NULL) const
    {out << "u_" << varIndex_;}

  /// return the value of the variable
  inline CouNumber operator () () 
    {return (currValue_ = expression::ubounds_ [varIndex_]);}

  /// differentiation
  inline expression *differentiate (int) 
    {return new exprConst (0);}

  /// dependence on variable set
  inline int dependsOn (int *, int) 
    {return 0;}

  /// get a measure of "how linear" the expression is:
  virtual inline int Linearity () 
    {return CONST;}

  /// code for comparisons
  virtual enum expr_type code ()
    {return COU_EXPRUBOUND;}
};

#endif
