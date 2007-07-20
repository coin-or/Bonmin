/*
 * Name:    exprLog.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of logarithm
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRLOG_H
#define COUENNE_EXPRLOG_H

#include <exprInv.hpp>
#include <expression.hpp>
#include <math.h>


// class logarithm

class exprLog: public exprUnary {

 public:

  // Constructors, destructor
  exprLog  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprLog (argument_ -> clone ());}

  /// the operator's function
  inline unary_function F () {return log;}

  /// print operator
  std::string printOp () const
    {return "log";}

  // differentiation
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // return expression of this same type with argument arg
  inline expression *mirror (expression *arg)
    {return new exprLog (arg);}

  // return derivative of univariate function of this type
  inline expression *mirror_d (expression *arg)
    {return new exprInv (arg);}

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRINV;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);
};

#endif
