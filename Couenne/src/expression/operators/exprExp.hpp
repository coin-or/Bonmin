/*
 * Name:    exprExp.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPREXP_H
#define COUENNE_EXPREXP_H

#include <exprUnary.hpp>
#include <math.h>

/// class for the exponential

class exprExp: public exprUnary {

 public:

  /// Constructor
  exprExp (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// Cloning method
  expression *clone () const
    {return new exprExp (argument_ -> clone ());}

  /// The operator's function
  inline unary_function F () {return exp;}

  /// Print operator
  std::string printOp () const
    {return "exp";}

  /// Differentiation
  expression *differentiate (int index); 

  /// Return expression of this same type with argument arg
  inline expression *mirror (expression *arg)
    {return new exprExp (arg);}

  /// Return derivative of univariate function of this type
  inline expression *mirror_d (expression *arg)
    {return new exprExp (arg);}

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Generate convexification cuts for this expression
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// Code for comparisons
  virtual enum expr_type code () {return COU_EXPREXP;}

  /// Implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);
};

#endif
