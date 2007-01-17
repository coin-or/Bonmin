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


// expression clone (points to VALUE and EXPRESSION of another expression) 

class exprClone: public exprCopy {

 public:

  // Constructor, destructor
  exprClone  (expression *copy): 
    exprCopy (copy -> Original ()) {}
  ~exprClone () {}

  // copy constructor
  exprClone (const exprClone &e):
    exprCopy (e.Original ()) {}

  // cloning method
  virtual exprClone *clone ()
  {return new exprClone (this);}

  // If this is an exprClone of a exprClone of an expr???, point to
  // the original expr??? instead of an exprClone -- improves computing
  // efficiency
  virtual inline expression *Original () const {return copy_ -> Original ();}

  // I/O
  void print (std::ostream &out) 
    {out << "{"; copy_ -> print (out); out << "}";}

  // value (empty)
  inline CouNumber Value () const 
    {return copy_ -> Value ();}

  // null function for evaluating the expression
  inline CouNumber operator () () 
    {return (currValue_ = (*copy_) ());}
};

#endif
