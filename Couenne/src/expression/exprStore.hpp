/*
 * Name:    exprStore.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of a storage class for expressions
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSTORE_HPP
#define COUENNE_EXPRSTORE_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "exprCopy.hpp"


/// storage class for previously evaluated expressions

class exprStore: public exprCopy {

 protected:

  /// Value of the (previously evaluated) expression
  CouNumber value_;

 public:

  /// Constructor
  exprStore (expression *copy):
    exprCopy (copy) {}

  /// Store constructor
  exprStore (const exprStore &e, const std::vector <exprVar *> *variables = NULL):
    exprCopy (e, variables) {
    //copy_  = e.Original () -> clone ();
  }

  /// Destructor
  virtual ~exprStore () 
  {copy_ = NULL;}

  /// Cloning method
  virtual expression *clone (const std::vector <exprVar *> *variables = NULL) const
  {return new exprStore (*this, variables);}

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
  {return (copy_ -> Value ());}
};

#endif
