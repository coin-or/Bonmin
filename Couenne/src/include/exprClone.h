/*
 * Name:    exprClone.h
 * Author:  Pietro Belotti
 * Purpose: definition of the clone class (different from exprCopy in
 *          that evaluation are propagated)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCLONE_H
#define COUENNE_EXPRCLONE_H

#include <iostream>

#include <CouenneTypes.h>
#include <exprCopy.h>


/// expression clone (points to VALUE and EXPRESSION of another expression) 

class exprClone: public exprCopy {

 public:

  /// Constructor
  exprClone  (expression *copy): 
    exprCopy (copy) {}

  /// destructor
  //  ~exprClone () {}

  /// copy constructor
  exprClone (const exprClone &e):
    exprCopy (e) {}

  /// cloning method
  exprClone *clone () const
  {return new exprClone (*this);}

  /// I/O
  void print (std::ostream &out) const
    //{out << "{"; copy_ -> Original () -> print (out); out << "}";}
    {copy_ -> Original () -> print (out);}
    //{copy_ -> print (out);}
    //{out << ",";}

  /// value
  inline CouNumber Value () const 
    {return copy_ -> Value ();}

  /// null function for evaluating the expression
  inline CouNumber operator () () 
    {return (currValue_ = (*copy_) ());}
};

#endif
