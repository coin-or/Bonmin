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
#include "exprVar.hpp"


/// variable-type operator. All variables of the expression must be
/// objects of this class

class exprIVar: public exprVar {

 public:

  /// Constructor
  exprIVar (int varIndex):
    exprVar (varIndex) {}

  /// Copy constructor
  exprIVar (const exprIVar &e, const std::vector <exprVar *> *variables = NULL):
    exprVar (e, variables) {}

  /// Cloning method
  virtual exprVar *clone (const std::vector <exprVar *> *variables = NULL) const
  {return ((variables && (*variables) [varIndex_]) ? 
	   (*variables) [varIndex_] :
	   new exprIVar (*this, variables));}

  //{return (variables ? (*variables) [varIndex_] : new exprIVar (*this, variables));}
  //{return (//keep_variables ? new exprClone (this) : 
  //new exprIVar (*this, variables));}

  /// Print
  virtual void print (std::ostream &out = std::cout, bool = false) const
    {out << "y_" << varIndex_;}

  /// Is this expression integer?
  virtual inline bool isInteger ()
    {return true;}
};

#endif
