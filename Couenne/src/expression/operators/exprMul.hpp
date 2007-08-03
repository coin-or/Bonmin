/*
 * Name:    exprMul.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of multiplications
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRMUL_H
#define COUENNE_EXPRMUL_H

#include <exprOp.hpp>
#include <exprAux.hpp>
#include <CouenneProblem.hpp>


/// class for multiplications

class exprMul: public exprOp {

 public:

  /// Constructor
  exprMul (expression **, int);

  /// constructor with two arguments
  exprMul (expression *, expression *);

  /// cloning method
  expression *clone () const
    {return new exprMul (clonearglist (), nargs_);}

  /// print operator
  std::string printOp () const
    {return "*";}

  /// method to evaluate the expression
  inline CouNumber operator () ();

  /// differentiation
  expression *differentiate (int index); 

  /// simplification
  expression *simplify ();

  /// get a measure of "how linear" the expression is:
  virtual int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  /// generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  expression *getFixVar ();

  /// code for comparison
  virtual enum expr_type code () {return COU_EXPRMUL;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);
};


/// compute multiplication

inline CouNumber exprMul:: operator () () {

  register CouNumber ret = 1;

  expression **al = arglist_;

  for (register int n = nargs_; n--;)
    ret *= (**al++) ();

  return (currValue_ = ret);
}


/// unified convexification of products and divisions

void unifiedProdCuts (const CouenneCutGenerator *, OsiCuts &, 
		      int, CouNumber, CouNumber, CouNumber,
		      int, CouNumber, CouNumber, CouNumber,
		      int, CouNumber, CouNumber, CouNumber,
		      t_chg_bounds *);

#endif
