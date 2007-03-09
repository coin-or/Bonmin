/*
 * Name:    exprIVar.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprIVar for integer variables 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRIVAR_H
#define COUENNE_EXPRIVAR_H

#include <iostream>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprConst.h>
#include <exprVar.h>

class CouenneProblem;


// variable-type operator. All variables of the expression must be
// objects of this class

class exprIVar: public exprVar {

 public:

  // Constructor
  exprIVar (int varIndex):
    exprVar (varIndex) {}

  // copy constructor
  exprIVar (const exprIVar &e):
    exprVar (e.Index ()) {}

  // cloning method
  virtual exprIVar *clone () const
    {return new exprIVar (*this);}

  // print
  virtual void print (std::ostream &out) const
    {out << "y_" << varIndex_;}

  // is this expression integer?
  virtual bool isInteger ()
    {return true;}
};

#endif
