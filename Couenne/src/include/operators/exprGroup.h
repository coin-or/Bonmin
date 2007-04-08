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

    if (index_) {
      delete [] index_;
      delete [] coeff_;
    }
  }

  /// get constant, indices and coefficients vectors
  CouNumber  getc0      () {return c0_;}
  CouNumber *getCoeffs  () {return coeff_;}
  int       *getIndices () {return index_;}

  /// cloning method
  virtual expression *clone () const
    {return new exprGroup (*this);}

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
  //  virtual exprAux *standardize (CouenneProblem *p);

  /// generate equality between *this and *w
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg);

  /// only compare with people of the same kind
  virtual int compare (exprGroup &);

  ///
  virtual enum expr_type code () {return COU_EXPRGROUP;}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ();
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprGroup::operator () () {

  register CouNumber  ret = c0_ + exprSum::operator () (),
                     *coe = coeff_;

  for (register int *ind = index_; *ind >= 0;)
    ret += *coe++ * expression::Variable (*ind++);

  return (currValue_ = ret);
}

#endif
