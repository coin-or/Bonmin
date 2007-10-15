/*
 * Name:    exprAux.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the auxiliary variable class (used in
 *          standardization and convexification)
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRAUX_HPP
#define COUENNE_EXPRAUX_HPP

#include <iostream>
#include <exprVar.hpp>
#include <exprBound.hpp>
#include <exprMax.hpp>
#include <exprMin.hpp>
#include <CouenneTypes.hpp>
#include <CglCutGenerator.hpp>


/** Auxiliary variable
 *
 *  It is associated with an expression which may depend on original
 *  and/or other auxiliary variables. It is used for AMPL's defined
 *  variables (aka common expressions) and to reformulate nonlinear
 *  constraints/objectives.
 *
 */

class exprAux: public exprVar {

 protected:

  /// The expression associated with this auxiliary variable
  expression *image_;

  /// lower bound, a function of the associated expression and the
  /// bounds on the variables in the expression
  expression *lb_;

  /// upper bound, a function of the associated expression and the
  /// bounds on the variables in the expression
  expression *ub_;

  /// used in rank-based branching variable choice: original variables
  /// have rank 1; auxiliary w=f(x) has rank r(w) = r(x)+1; finally,
  /// auxiliary w=f(x1,x2...,xk) has rank r(w) = 1+max{r(xi):i=1..k}.
  int rank_;

  /// number of appearances of this aux in the formulation. The more
  /// times it occurs in the formulation, the more implication its
  /// branching has on other variables
  int multiplicity_;

  /// is this variable integer?
  bool integer_;

 public:

  /// Node type
  inline enum nodeType Type () 
    {return AUX;}

  /// Constructor
  exprAux (expression *, int, int, bool = false);

  /// Constructor to be used with standardize ([...], false)
  exprAux (expression *);

  /// Destructor
  ~exprAux ();

  /// Copy constructor
  exprAux (const exprAux &);

  /// Cloning method
  virtual exprAux *clone () const
    {return new exprAux (*this);}

  expression *Lb () {return lb_;} ///< get lower bound expression
  expression *Ub () {return ub_;} ///< get upper bound expression

  /// Print expression
  virtual void print (std::ostream & = std::cout, 
		      bool = false, CouenneProblem * = NULL) const;

  /// The expression associated with this auxiliary variable
  inline expression *Image () const
    {return image_;}

  /// Null function for evaluating the expression
  inline CouNumber operator () () 
    {return (currValue_ = expression::Variable (varIndex_));}

  /// Differentiation
  inline expression *differentiate (int index) 
    {return image_ -> differentiate (index);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  int DepList (std::set <int> &deplist, 
	       enum dig_type type = ORIG_ONLY,
	       CouenneProblem *p = NULL);

  /// simplify
  expression *simplify ();

  /// Get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return LINEAR;}
    /*return image_ -> Linearity ();*/

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&lb, expression *&ub);

  /// set bounds depending on both branching rules and propagated
  /// bounds. To be used after standardization
  void crossBounds ();

  /// generate cuts for expression associated with this auxiliary
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, const CouenneCutGenerator *, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// used in rank-based branching variable choice
  virtual inline int rank (CouenneProblem *p = NULL)
    {return rank_;} 

  /// is this expression integer?
  virtual inline bool isInteger () {
    return ((integer_ == AUX_INTEGER) || 
	    (integer_ == AUX_UNSET) && ((integer_ = image_ -> isInteger ()) == AUX_INTEGER));
  }

  /// Tell this variable appears once more
  inline void increaseMult () {++multiplicity_;}

  /// Tell this variable appears once less (standardized within
  /// exprSum, for instance)
  inline void decreaseMult () {--multiplicity_;}

  /// How many times this variable appears 
  inline int Multiplicity () {return multiplicity_;}
};

/// allow to draw function within intervals and cuts introduced
void draw_cuts (OsiCuts &, const CouenneCutGenerator *, 
		int, expression *, expression *);

#endif
