/*
 * Name:    exprInv.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of inverse of a function (1/f(x))
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRINV_H
#define COUENNE_EXPRINV_H

#include <exprUnary.hpp>


/// the operator itself
inline CouNumber inv (register CouNumber arg) 
{return 1.0 / arg;}


/// derivative of inv (x)
inline CouNumber oppInvSqr (register CouNumber x) 
{return (- inv (x*x));}


/// inv_dblprime, second derivative of inv (x)
inline CouNumber inv_dblprime (register CouNumber x) 
{return (2 * inv (x*x*x));}


/// class inverse (1/f(x))

class exprInv: public exprUnary {

 public:

  /// Constructors, destructor
  exprInv  (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone () const
    {return new exprInv (argument_ -> clone ());}

  /// the operator's function
  inline unary_function F () {return inv;}

  /// output "1/argument"
  virtual void print (std::ostream &out = std::cout, bool = false, CouenneProblem * = NULL) const;

  /// differentiation
  expression *differentiate (int index); 

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () {
    if (argument_ -> Type () == CONST) return CONSTANT;
    else                               return NONLINEAR;
  }

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// generate equality between *this and *w
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
