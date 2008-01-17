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

#include <vector>

#include "exprSum.hpp"
#include "exprVar.hpp"


/// class Group, with constant, linear and nonlinear terms

class exprGroup: public exprSum {

public:

  typedef std::vector <std::pair <exprVar *, CouNumber> > lincoeff;

protected:

  mutable lincoeff lcoeff_;  ///< coefficients and indices of the linear term
  CouNumber        c0_;      ///< constant term

public:

  /// Constructor
  exprGroup  (CouNumber,
	      lincoeff &, 
	      expression ** = NULL, 
	      int = 0);

  /// Copy constructor
  exprGroup (const exprGroup &src);

  // Get constant, indices, and coefficients vectors, and number of linear terms
  CouNumber  getc0      () {return c0_;}      ///< return constant term

  lincoeff &lcoeff () const {return lcoeff_;} ///< return linear term coefficients

  /// Cloning method
  virtual expression *clone () const
    {return new exprGroup (*this);}

  /// Print expression to iostream
  virtual void print (std::ostream &   = std::cout, 
		      bool             = false) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual int DepList (std::set <int> &deplist, 
		       enum dig_type type = ORIG_ONLY);

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
  virtual void generateCuts (expression *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// only compare with people of the same kind
  virtual int compare (exprGroup &);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRGROUP;}

  /// is this expression integer?
  virtual bool isInteger ();

  /// used in rank-based branching variable choice
  virtual int rank ();

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);

  /// replace variable x with new (aux) w
  virtual void replace (exprVar *x, exprVar *w);
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprGroup::operator () () {

  CouNumber ret = c0_ + exprSum::operator () (); // add constant and nonlinear part

  // add linear part
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    ret += el -> second * (*(el -> first)) ();

  return (CouNumber) ret;
}

#endif
