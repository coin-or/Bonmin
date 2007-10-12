/*
 * Name:    exprPow.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRPOW_HPP
#define COUENNE_EXPRPOW_HPP

#include <math.h>

#include <exprOp.hpp>
#include <exprMul.hpp>
#include <exprClone.hpp>
#include <exprConst.hpp>


/// Power of an expression (binary operator)

class exprPow: public exprOp {

 public:

  /// Constructor
  exprPow  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// Constructor with only two arguments
  exprPow (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  /// cloning method
  expression *clone () const
    {return new exprPow (clonearglist (), nargs_);}

  /// print operator
  std::string printOp () const
    {return "^";}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// differentiation
  expression *differentiate (int index); 

  /// simplification
  expression *simplify ();

  /// get a measure of "how linear" the expression is
  int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  exprAux *standardize (CouenneProblem *p, bool addAux = true);

  /// generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  expression *getFixVar () 
    {return arglist_ [0];}

  /// code for comparison
  virtual enum expr_type code () 
    {return COU_EXPRPOW;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);
};


/// compute power and check for integer-and-odd inverse exponent

inline CouNumber safe_pow (register CouNumber base, 
			   register CouNumber exponent) {

  if (base < 0) {

    register int rndexp;

    if (((fabs (exponent - (rndexp = COUENNE_round (exponent))) < COUENNE_EPS) ||
	 ((fabs (exponent) > COUENNE_EPS) && 
	  (fabs (1. / exponent - (rndexp = COUENNE_round (1. / exponent))) < COUENNE_EPS)))
	&& (rndexp % 2))
      return (- pow (- base, exponent));
  }

  if (fabs (base) >= COUENNE_INFINITY) {

    if (base <= -COUENNE_INFINITY) {

      register int intk = COUENNE_round (exponent);

      if ((fabs (exponent - intk) < COUENNE_EPS) && (intk % 2))
	return (exponent < 0) ? 0 : -COUENNE_INFINITY;
    }
    else return (exponent < 0) ? 0 : COUENNE_INFINITY;
  }

  return (pow (base, exponent));
}


/// compute power
inline CouNumber exprPow::operator () () {

  //  exprOp:: operator () ();

  //  register CouNumber exponent = *sp--;
  //  register CouNumber base     = *sp--;

  //  return (currValue_ = safe_pow (base, exponent));
  return (currValue_ = safe_pow ((**arglist_) (), (*(arglist_ [1])) ()));
}


/// add upper/lower envelope to power in convex/concave areas
void addPowEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int,
		     CouNumber, CouNumber, CouNumber, CouNumber, CouNumber, int);


/// find proper tangent point to add deepest tangent cut
CouNumber powNewton (CouNumber, CouNumber, unary_function, unary_function, unary_function);

#endif
