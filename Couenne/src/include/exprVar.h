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

  int varIndex_; //< the address of the variable's current value

 public:

  // node type
  virtual inline enum nodeType Type () 
    {return VAR;}

  // Constructor, destructor
  exprVar  (int varIndex):
    varIndex_ (varIndex) {}

  ~exprVar () {}

  // get variable index in problem
  inline int Index () 
    {return varIndex_;}

  // print
  virtual void print (std::ostream &out) 
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

  // Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *w, int *&, expression ***&, int **&, 
  //		       expression **&, enum con_sign *&);

  // similarly, construct linear over-estimator for expression within
  // problem *p (p is used to add convexification constraints). It is
  // also used when this function appears with a minus sign in the
  // expression
  //  inline int upperLinearHull (exprAux *, int *&, expression ***&, int **&, 
  //		       expression **&, enum con_sign *&)
  //    {return 0;}

  // generate convexification cut for constraint w = this
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
