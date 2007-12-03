/*
 * Name:    exprGroup.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of mixed sum expressions (constant+linear+nonlinear)
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRGROUP_H
#define COUENNE_EXPRGROUP_H

#include <exprSum.hpp>


/// class Group, with constant, linear and nonlinear terms

class exprGroup: public exprSum {

 protected:

  int        nlterms_; ///< number of linear terms
  int       *index_;   ///< indices of linear terms (terminated by a -1)
  CouNumber *coeff_;   ///< coefficient of linear terms
  CouNumber  c0_;      ///< constant term

 public:

  /// Constructor
  exprGroup  (CouNumber, int *, CouNumber *, expression ** = NULL, int = 0);

  /// Copy constructor
  exprGroup (const exprGroup &src);

  /// Destructor
  virtual ~exprGroup ();

  // Get constant, indices, and coefficients vectors, and number of linear terms
  CouNumber  getc0      () {return c0_;}      ///< return constant term
  CouNumber *getCoeffs  () {return coeff_;}   ///< return coefficients of linear part
  int       *getIndices () {return index_;}   ///< return indices of linear part
  int        getnLTerms () {return nlterms_;} ///< return number of linear terms

  /// Cloning method
  virtual expression *clone () const
    {return new exprGroup (*this);}

  /// Print expression to iostream
  virtual void print (std::ostream & = std::cout, bool = false, 
		      CouenneProblem * = NULL) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual int DepList (std::set <int> &deplist, 
		       enum dig_type type = ORIG_ONLY,
		       CouenneProblem *p = NULL);

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

  /// is this expression integer?  TODO: find a way to get integrality
  /// of linear terms, as for now we have to return false
  virtual inline bool isInteger ()
    {return false;}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprGroup::operator () () {

  register CouNumber
     ret  = c0_ + exprSum::operator () (), // add constant and nonlinear part
    *coe  = coeff_,
    *vars = expression::Variables ();

  // add linear part
  for (register int n = nlterms_, *ind = index_; n--;)
    ret += *coe++ * vars [*ind++];

  return ret;
}


/// scan subexpression for single index
int scanIndex (int index, std::set <int> &deplist, CouenneProblem *p, enum dig_type type);

#endif
