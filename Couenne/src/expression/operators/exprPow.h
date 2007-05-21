/*
 * Name:    exprPow.h
 * Author:  Pietro Belotti
 * Purpose: definition of powers
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRPOW_H
#define COUENNE_EXPRPOW_H

#include <math.h>

#include <exprOp.h>
#include <exprMul.h>
#include <exprClone.h>
#include <exprConst.h>


// Power of an expression (binary operator)

class exprPow: public exprOp {

 public:

  // Constructors, destructor
  exprPow  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprPow (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  // cloning method
  expression *clone () const
    {return new exprPow (clonearglist (), nargs_);}

  /// print operator
  std::string printOp () const
    {return "^";}

  // I/O
  //  void print (std::ostream&) const;

  // function for the evaluation of the expression
  CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is
  int Linearity ();

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  expression *getFixVar () 
    {return arglist_ [0];}

  /// code for comparison
  virtual enum expr_type code () 
    {return COU_EXPRPOW;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, char *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);

  /*  /// distance covered by current point if branching rule applied to this expression
  double BranchGain (expression *, const OsiBranchingInformation *);

  /// branching object best suited for this expression
  OsiBranchingObject *BranchObject (expression *, const OsiBranchingInformation *);*/
};


// compute power and check for integer-and-odd inverse exponent

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


// compute power
inline CouNumber exprPow::operator () () {

  exprOp:: operator () ();

  register CouNumber exponent = *sp--;
  register CouNumber base     = *sp--;

  return (currValue_ = safe_pow (base, exponent));
}


// add upper/lower envelope to power in convex/concave areas
void addPowEnvelope (const CouenneCutGenerator *, OsiCuts &, int, int,
		     CouNumber, CouNumber, CouNumber, CouNumber, CouNumber, int);


// find proper tangent point to add deepest tangent cut
CouNumber powNewton (CouNumber, CouNumber, unary_function, unary_function, unary_function);

#endif
