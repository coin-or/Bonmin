/*
 * Name:    exprSub.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of subtractions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUB_H
#define COUENNE_EXPRSUB_H

#include <exprOp.hpp>
#include <CouennePrecisions.h>
#include <CouenneProblem.hpp>

/// class for subtraction

class exprSub: public exprOp {

 public:

  /// Constructors, destructor
  exprSub  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprSub (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  /// cloning method
  expression *clone () const
    {return new exprSub (clonearglist (), nargs_);}

  //// print operator
  std::string printOp () const
    {return "-";}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// differentiation
  expression *differentiate (int index); 

  /// simplification
  expression *simplify ();

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () {

    int lin1 = arglist_ [0] -> Linearity ();
    int lin2 = arglist_ [1] -> Linearity ();

    if (lin1 < lin2) return lin2;
    else             return lin1;
  }

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  /// special version for linear constraints
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1,
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRSUB;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);
};


/// compute subtraction

inline CouNumber exprSub::operator () ()
{return (currValue_ = (*(*arglist_)) () - (*(arglist_ [1])) ());}

#endif
