/*
 * Name:    exprAbs.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the absolute value of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRABS_H
#define COUENNE_EXPRABS_H

#include <math.h>

#include <exprUnary.hpp>
#include <exprConst.hpp>


/// class for \f$w=|f(x)|\f$

class exprAbs: public exprUnary {

 public:

  /// Constructor
  exprAbs  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// The operator's function
  inline unary_function F () {return fabs;}

  /// cloning method
  expression *clone () const
    {return new exprAbs (argument_ -> clone ());}

  /// output
  std::string printOp () const
    {return "abs";}

  /// differentiation
  expression *differentiate (int index); 

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  enum expr_type code () {return COU_EXPRABS;}

  /// is this expression integer?
  bool isInteger ()
  {return argument_ -> isInteger ();}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);
};

#endif
