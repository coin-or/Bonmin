/*
 * Name:    exprVar.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprVar for variables 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRVAR_H
#define COUENNE_EXPRVAR_H

#include <iostream>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprConst.h>

class CouenneProblem;


// variable-type operator. All variables of the expression must be
// objects of this class

class exprVar: public expression {

 protected:

  int varIndex_; //< the index of the variable's current value

 public:

  // node type
  virtual inline enum nodeType Type () 
    {return VAR;}

  // Constructor
  exprVar (int varIndex):
    varIndex_ (varIndex) {}

  // destructor
  virtual  ~exprVar () {}

  // copy constructor
  exprVar (const exprVar &e):
    varIndex_ (e.Index ()) {}

  // cloning method
  virtual exprVar *clone ()
    {return new exprVar (*this);}

  // get variable index in problem
  inline int Index () const
    {return varIndex_;}

  // string equivalent
  virtual const std::string name () const;

  // print
  virtual void print (std::ostream &out) const
    {out << "x_" << varIndex_;}

  // return the value of the variable
  virtual inline CouNumber operator () () 
    {return (currValue_ = expression::variables_ [varIndex_]);}

  // return the value of the variable
  inline CouNumber Value ()
    {return currValue_;}

  // differentiation
  virtual inline expression *differentiate (int index) 
    {return new exprConst ((index == varIndex_) ? 1 : 0);}

  // dependence on variable set
  virtual bool dependsOn (int *, int);

  // simplify
  inline expression *simplify () 
    {return NULL;}

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  virtual inline int Linearity ()
    {return LINEAR;}

  // is this expression integer?
  virtual bool isInteger ()
    {return false;}

  // Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  // generate convexification cut for constraint w = this
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar () {return this;}
};

#endif
