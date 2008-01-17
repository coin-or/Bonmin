/*
 * Name:    exprConst.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprConst
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCONST_HPP
#define COUENNE_EXPRCONST_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprClone.hpp"


/// constant-type operator

class exprConst: public expression {

private: 

  /// the value of this constant
  CouNumber value_;

public:

  /// node type
  inline enum nodeType Type () 
    {return CONST;}

  /// value of expression
  inline CouNumber Value () const 
    {return value_;}

  /// Constructor
  exprConst (CouNumber value)
    {value_ = value;}

  /// Copy constructor
  exprConst (const exprConst &e)
    {value_ = e.value_;}

  /// Cloning method
  virtual exprConst *clone () const
    {return new exprConst (value_);}

  /// I/O
  void print (std::ostream &out = std::cout, 
	      bool = false) const
    {out << value_;}

  /// return constant's value
  inline CouNumber operator() () 
    {return value_;}

  /// differentiation
  inline expression *differentiate (int) 
    {return new exprConst (0.);}

  /// dependence on variable set
  int dependsOn (int *ind, int n, enum dig_type type = STOP_AT_AUX)
    {return 0;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return ((fabs (value_) < COUENNE_EPS) ? ZERO: CONSTANT);}

  /// Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) {
    lower = new exprClone (this);
    upper = new exprClone (this);
  }

  /// generate convexification cut for constraint w = this
  void generateCuts (expression *, const OsiSolverInterface &, 
		     OsiCuts &, const CouenneCutGenerator *, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code () 
  {return COU_EXPRCONST;}

  /// is this expression integer?
  virtual bool isInteger () 
  {return ::isInteger (value_);}

  /// used in rank-based branching variable choice
  virtual int rank ()
    {return 0;} 
};

#endif
