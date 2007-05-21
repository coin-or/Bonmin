/*
 * Name:    exprExp.h
 * Author:  Pietro Belotti
 * Purpose: definition of the exponential of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPREXP_H
#define COUENNE_EXPREXP_H

#include <exprUnary.h>
#include <math.h>

// class for the exponential

class exprExp: public exprUnary {

 public:

  // Constructors, destructor
  exprExp  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprExp (argument_ -> clone ());}

  /// the operator's function
  inline unary_function F () {return exp;}

  /// print operator
  std::string printOp () const
    {return "exp";}

  /// print position
  //  virtual const std::string printPos () 
  //    {return PRE;}


  // output
  //  void print (std::ostream&) const;

  // differentiation
  expression *differentiate (int index); 

  // return expression of this same type with argument arg
  inline expression *mirror (expression *arg)
    {return new exprExp (arg);}

  // return derivative of univariate function of this type
  inline expression *mirror_d (expression *arg)
    {return new exprExp (arg);}

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPREXP;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, char *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);

  /*/// distance covered by current point if branching rule applied to this expression
  double BranchGain (expression *, const OsiBranchingInformation *);

  /// branching object best suited for this expression
  OsiBranchingObject *BranchObject (expression *, const OsiBranchingInformation *);*/
};

#endif
