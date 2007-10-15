/*
 * Name:    exprSub.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of subtractions
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUB_HPP
#define COUENNE_EXPRSUB_HPP

#include <exprOp.hpp>
#include <CouennePrecisions.hpp>
#include <CouenneProblem.hpp>

/// class for subtraction

class exprSub: public exprOp {

 public:

  /// Constructor
  exprSub  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// Constructor with two explicit elements
  exprSub (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  /// Cloning method
  expression *clone () const
    {return new exprSub (clonearglist (), nargs_);}

  //// Print operator
  std::string printOp () const
    {return "-";}

  /// Function for the evaluation of the difference
  CouNumber operator () ();

  /// Differentiation
  expression *differentiate (int index); 

  /// Simplification
  expression *simplify ();

  /// Get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () {

    int lin1 = arglist_ [0] -> Linearity ();
    int lin2 = arglist_ [1] -> Linearity ();

    if (lin1 < lin2) return lin2;
    else             return lin1;
  }

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// Reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// Special version for linear constraints
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// is this expression integer?
  virtual bool isInteger ();

  /// Code for comparisons
  virtual enum expr_type code () {return COU_EXPRSUB;}

  /// Implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);
};


/// Compute difference

inline CouNumber exprSub::operator () ()
{return (currValue_ = (*(*arglist_)) () - (*(arglist_ [1])) ());}

#endif
