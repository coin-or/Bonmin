/*
 * Name:    exprSum.h
 * Author:  Pietro Belotti
 * Purpose: definition of sum expressions
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSUM_H
#define COUENNE_EXPRSUM_H

#include <exprOp.h>


///  class sum 

class exprSum: public exprOp {

 public:

  /// Constructors, destructor
  exprSum  (expression **, int);

  /// Constructor with two elements
  exprSum (expression *, expression *);

  /// cloning method
  virtual expression *clone () const
    {return new exprSum (clonearglist (), nargs_);}

  /// I/O
  virtual void print (std::ostream &) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// differentiation
  virtual expression *differentiate (int index); 

  /// simplification
  virtual expression *simplify ();

  /// get a measure of "how linear" the expression is:
  virtual int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  virtual exprAux *standardize (CouenneProblem *p);

  /// generate equality between *this and *w
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg);

  ///
  virtual enum expr_type code () {return COU_EXPRSUM;}

  /// implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, char *);

  ///
  virtual expression *getFixVar () {return *arglist_;}
};


/// compute sum

inline CouNumber exprSum::operator () () {

  exprOp:: operator () ();

  register CouNumber ret = *sp--; 
  register int       n   = nargs_;

  while (--n)
    ret += *sp--;

  return (currValue_ = ret);
}

#endif
