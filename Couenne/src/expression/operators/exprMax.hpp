/*
 * Name:    exprMax.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of $\f(x_{\argmax_{i\in I} y_i})$ 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRMAX_H
#define COUENNE_EXPRMAX_H

#include "exprOp.hpp"
#include "exprCopy.hpp"
#include "exprStore.hpp"


/// class for maxima

class exprMax: public exprOp {

 public:

  /// Constructor
  exprMax  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// Constructor with only two arguments
  exprMax  (expression *el0, expression *el1):
    exprOp (new expression * [4], 4) {
    arglist_ [0] = new exprCopy (el0); arglist_ [1] = new exprStore (arglist_ [0]);
    arglist_ [2] = new exprCopy (el1); arglist_ [3] = new exprStore (arglist_ [2]);
  }

  /// cloning method
  exprMax *clone () const
    {return new exprMax (clonearglist (), nargs_);}

  /// print operator
  std::string printOp () const
    {return "max";}

  /// print position
  enum pos printPos () const
    {return PRE;}

  /// function for the evaluation of the expression
  inline CouNumber operator () ();

  /// differentiation
  inline expression *differentiate (int) 
    {return NULL;} 

  /// simplification
  inline expression *simplify () 
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () 
    {return NONLINEAR;}

  // Get lower and upper bound of an expression (if any)
  //  void getBounds (expression *&, expression *&);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual inline exprAux *standardize (CouenneProblem *, bool addAux = true)
    {return NULL;}

  /// generate equality between *this and *w
  void generateCuts (expression *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code ()
  {return COU_EXPRMAX;}
};


/// compute maximum

inline CouNumber exprMax::operator () () {

  CouNumber best_val = (*(arglist_ [0])) ();
  int best_ind = 0;

  for (register int ind = 2; ind < nargs_; ind += 2) {

    register CouNumber val = (*(arglist_ [ind])) ();

    if (val > best_val) {
      best_ind = ind;
      best_val = val;
    }
  }

  return (*(arglist_ [best_ind + 1])) ();
}

#endif
