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

  /// Store constructor -- Must go
  exprStore (const exprStore &e, Domain *d = NULL):
    exprCopy (e, d) {
    //copy_  = e.Original () -> clone ();
  }

  /// Destructor
  virtual ~exprStore () 
  {copy_ = NULL;}

  /// I/O -- Must go
  virtual void print (std::ostream &out = std::cout, 
		      bool descend      = false) const
  {out << "<"; copy_ -> print (out, descend); out << ">"; }

  /// Cloning method
  virtual inline expression *clone (Domain *d = NULL) const
  {return new exprStore (*this, d);}

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
  {return (copy_ -> Value ());}

  /// empty function to redirect variables to proper variable vector
  virtual void realign (const CouenneProblem *p);
};

#endif
