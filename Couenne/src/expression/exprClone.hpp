/*
 * Name:    exprClone.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the clone class (different from exprCopy in
 *          that evaluation are propagated)
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCLONE_HPP
#define COUENNE_EXPRCLONE_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "exprCopy.hpp"


/// expression clone (points to another expression) 

class exprClone: public exprCopy {

 public:

  /// Constructor
  exprClone  (expression *copy): 
    exprCopy (copy) {}

  /// copy constructor
  exprClone (const exprClone &e, Domain *d = NULL):
    exprCopy (e, d) {}

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprClone (*this, d);}

  /// Destructor
  virtual ~exprClone () 
  {copy_ = NULL;}

  /// I/O
  //void print (std::ostream &out = std::cout, 
  //bool descend      = false) const
    //{copy_ -> Original () -> print (out, descend);}
  //{out << "{"; copy_ -> print (out, descend); out << "}"; } // Must go

  /// value
  inline CouNumber Value () const 
  {return copy_ -> Value ();}

  /// null function for evaluating the expression
  inline CouNumber operator () () 
  {return ((*copy_) ());}

  /// empty function to redirect variables to proper variable vector
  void realign (const CouenneProblem *p);
};

#endif
