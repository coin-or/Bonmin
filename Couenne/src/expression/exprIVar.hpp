/*
 * Name:    exprIVar.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprIVar for integer variables 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRIVAR_HPP
#define COUENNE_EXPRIVAR_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprConst.hpp"
#include "exprVar.hpp"

class CouenneProblem;


/// variable-type operator. All variables of the expression must be
/// objects of this class

class exprIVar: public exprVar {

 public:

  /// Constructor
  exprIVar (int varIndex):
    exprVar (varIndex) {}

  /// copy constructor
  exprIVar (const exprIVar &e):
    exprVar (e.Index ()) {}

  /// cloning method
  virtual exprIVar *clone () const
    {return new exprIVar (*this);}

  /// print
  virtual void print (std::ostream &out = std::cout, bool = false, 
		      CouenneProblem * = NULL) const
    {out << "y_" << varIndex_;}

  /// is this expression integer?
  virtual inline bool isInteger ()
    {return true;}
};

#endif
