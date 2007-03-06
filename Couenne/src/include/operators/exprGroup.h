/*
 * Name:    exprGroup.h
 * Author:  Pietro Belotti
 * Purpose: definition of mixed sum expressions (constant+linear+nonlinear)
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRGROUP_H
#define COUENNE_EXPRGROUP_H

#include <exprSum.h>


///  class Group, with constant, linear and nonlinear terms

class exprGroup: public exprSum {

 protected:

  CouNumber  c0_;    //< constant term
  int       *index_; //< indices of linear terms (terminated by a -1)
  CouNumber *coeff_; //< coefficient of linear terms

 public:

  /// Constructor
  exprGroup  (CouNumber, int *, CouNumber *, expression ** = NULL, int = 0);

  /// copy constructor
  exprGroup (const exprGroup &src);

  /// destructor
  virtual ~exprGroup () {

    delete index_;
    delete coeff_;
  }

  /// cloning method
  virtual expression *clone () const
    {return new exprGroup (*this);}

  /// String equivalent (for comparisons)
  virtual const std::string name () const;

  /// I/O
  virtual void print (std::ostream &) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// differentiation
  virtual expression *differentiate (int index); 

  /// simplification
  virtual expression *simplify ()
    {exprOp::simplify (); return NULL;}

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
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprGroup::operator () () {

  register CouNumber ret = c0_ + exprSum::operator () ();

  for (register int *ind = index_, i=0; *ind >= 0;)
    ret += coeff_ [i++] * expression::Variable (*ind++);

  return (currValue_ = ret);
}

#endif
