/*
 * Name:    exprSub.h
 * Author:  Pietro Belotti
 * Purpose: definition of subtractions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUB_H
#define COUENNE_EXPRSUB_H

#include <exprOp.h>
#include <CouennePrecisions.h>
#include <CouenneProblem.h>

// class for subtraction

class exprSub: public exprOp {

 public:

  // Constructors, destructor
  exprSub  (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprSub (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  // cloning method
  expression *clone () const
    {return new exprSub (clonearglist (), nargs_);}

  // I/O
  void print (std::ostream&) const;

  // function for the evaluation of the expression
  CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () {

    int lin1 = arglist_ [0] -> Linearity ();
    int lin2 = arglist_ [1] -> Linearity ();

    if (lin1 < lin2) return lin2;
    else             return lin1;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  ///
  virtual enum expr_type code () {return COU_EXPRSUB;}
};


// compute subtraction

inline CouNumber exprSub::operator () () {

  exprOp:: operator () ();

  register CouNumber ret = - *sp--;
  return (currValue_ = ret + *sp--);
}

#endif
