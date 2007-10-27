/*
 * Name:    exprSum.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUM_H
#define COUENNE_EXPRSUM_H

#include "exprOp.hpp"


/// class sum 

class exprSum: public exprOp {

 public:

  /// Constructors, destructor
  exprSum  (expression **, int);

  /// Constructor with two elements
  exprSum (expression *, expression *);

  /// Empty destructor
  ~exprSum () {}
 
  /// Cloning method
  virtual expression *clone () const
    {return new exprSum (clonearglist (), nargs_);}

  /// Print operator
  std::string printOp () const
    {return "+";}

  /// Function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// Differentiation
  virtual expression *differentiate (int index); 

  /// Simplification
  virtual expression *simplify ();

  /// Get a measure of "how linear" the expression is:
  virtual int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// Special version for linear constraints
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// Code for comparison
  virtual enum expr_type code () 
    {return COU_EXPRSUM;}

  /// Implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// Checks for quadratic terms in the expression and returns an
  /// exprQuad if there are enough to create something that can be
  /// convexified
  exprAux *createQuadratic (CouenneProblem *);
};


/// compute sum

inline CouNumber exprSum::operator () () {

  register CouNumber ret = 0;

  expression **al = arglist_;

  for (register int n = nargs_; n--;)
    ret += (**al++) ();

  return (currValue_ = ret);
}

#endif
