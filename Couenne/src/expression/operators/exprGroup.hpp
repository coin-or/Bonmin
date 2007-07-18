/*
 * Name:    exprGroup.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of mixed sum expressions (constant+linear+nonlinear)
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRGROUP_H
#define COUENNE_EXPRGROUP_H

#include <exprSum.hpp>


///  class Group, with constant, linear and nonlinear terms

class exprGroup: public exprSum {

 protected:

  CouNumber  c0_;    //< constant term
  int       *index_; //< indices of linear terms (terminated by a -1)
  CouNumber *coeff_; //< coefficient of linear terms

  int nlterms_;      //< number of linear terms

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

  /// get constant, indices, and coefficients vectors, and number of linear terms
  CouNumber  getc0      () {return c0_;}
  CouNumber *getCoeffs  () {return coeff_;}
  int       *getIndices () {return index_;}
  int        getNLTerms () {return nlterms_;}

  /// cloning method
  virtual expression *clone () const
    {return new exprGroup (*this);}

  /// I/O
  virtual void print (std::ostream & = std::cout, bool = false, CouenneProblem * = NULL) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// check if exprGroup depends on a list of variables specified as parameters
  virtual int dependsOn (int * = NULL, int = 1);

  /// specialized version that checks variables in linear term through
  /// their image (if any) within the CouenneProblem
  virtual int dependsOn (CouenneProblem *, int * = NULL, int = 1);

  /// differentiation
  virtual expression *differentiate (int index); 

  /// simplification
  virtual expression *simplify ()
    {exprOp::simplify (); return NULL;}

  /// get a measure of "how linear" the expression is:
  virtual int Linearity ();

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// special version for linear constraints
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// only compare with people of the same kind
  virtual int compare (exprGroup &);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRGROUP;}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprGroup::operator () () {

  register CouNumber
     ret  = c0_ + exprSum::operator () (),
    *coe  = coeff_,
    *vars = expression::Variables ();

  for (register int *ind = index_; *ind >= 0;)
    ret += *coe++ * vars [*ind++];

  return (currValue_ = ret);
}

#endif
