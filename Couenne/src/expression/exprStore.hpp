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
  exprStore (const exprStore &e):
    exprCopy (e) {
    //copy_  = e.Original () -> clone ();
  }

  /// Cloning method
  virtual exprStore *clone () const
    {return new exprStore (*this);}

  /// value (the saved one)
  //  virtual inline CouNumber Value () const 
  //    {return value_;}

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
    {return (copy_ -> Value ());}
};

#endif
