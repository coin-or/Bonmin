/*
 * Name:    exprIVar.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprIVar for integer variables 
 *
 * (C) Carnegie-Mellon University, 2006-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRIVAR_HPP
#define COUENNE_EXPRIVAR_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprVar.hpp"


/// variable-type operator. All variables of the expression must be
/// objects of this class

class exprIVar: public exprVar {

 public:

  /// Constructor
  exprIVar (int varIndex, Domain *d = NULL):
    exprVar (varIndex, d) {}

  /// Copy constructor -- must go
  exprIVar (const exprIVar &e, Domain *d = NULL):
    exprVar (e, d) {}

  /// Cloning method
  virtual exprVar *clone (Domain *d = NULL) const
  {return new exprIVar (*this, d);}

  /// Print
  virtual void print (std::ostream &out = std::cout, bool = false) const
  {out << "y_" << varIndex_;}

  /// is this expression defined as an integer?
  virtual inline bool isDefinedInteger ()
  {return true;}

  /// Is this expression integer?
  virtual inline bool isInteger ()
  {return true;}
};

#endif
