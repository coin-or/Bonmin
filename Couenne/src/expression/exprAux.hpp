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

#include "expression.hpp"
#include "exprVar.hpp"

class CouenneCutGenerator;

/** Auxiliary variable
 *
 *  It is associated with an expression which may depend on original
 *  and/or other auxiliary variables. It is used for AMPL's defined
 *  variables (aka common expressions) and to reformulate nonlinear
 *  constraints/objectives.
 *
 */

class exprAux: public exprVar {

 public:

  /// integrality type of an auxiliary variable: unset, continuous, integer
  enum intType {Unset=-1, Continuous, Integer};

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
  enum intType integer_;

 public:

  /// Node type
  inline enum nodeType Type () const
    {return AUX;}

  /// Constructor
  exprAux (expression *, int, int, intType = Unset, Domain * = NULL);

  /// Constructor to be used with standardize ([...], false)
  exprAux (expression *, Domain * = NULL);

  /// Destructor
  ~exprAux ();

  /// Copy constructor
  exprAux (const exprAux &, Domain *d = NULL);

  /// Cloning method
  virtual inline exprVar *clone (Domain *d = NULL) const
  {return new exprAux (*this, d);}

  inline expression *Lb () {return lb_;} ///< get lower bound expression
  inline expression *Ub () {return ub_;} ///< get upper bound expression

  /// Print expression
  virtual void print (std::ostream & = std::cout, 
		      bool = false) const;

  /// The expression associated with this auxiliary variable
  inline expression *Image () const
    {return image_;}

  /// The expression associated with this auxiliary variable
  inline void Image (expression *image)
    {image_ = image;}

  /// Null function for evaluating the expression
  inline CouNumber operator () ()
    {return domain_ -> x (varIndex_);}

  /// Differentiation
  inline expression *differentiate (int index) 
    {return image_ -> differentiate (index);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  int DepList (std::set <int> &deplist, 
	       enum dig_type type = ORIG_ONLY);

  /// simplify
  expression *simplify ();

  /// Get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return LINEAR;}
    /*return image_ -> Linearity ();*/

  /// Get lower and upper bound of an expression (if any)
  //virtual void getBounds (expression *&lb, expression *&ub);

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
  virtual inline int rank ()
    {return rank_;} 

  /// is this expression defined as integer?
  virtual inline bool isDefinedInteger () {

    return ((integer_ == Integer) || 
	    (integer_ == Unset) && 
	    ((integer_ = (image_ -> isInteger ()) ? 
	      Integer : Continuous) == Integer));
  }

  /// is this expression integer?
  virtual inline bool isInteger () {

    if (isDefinedInteger ()) 
      return true;

    CouNumber l = lb ();
    return (::isInteger (l) && (fabs (l - ub ()) < COUENNE_EPS));
    //CouNumber l = (*(Lb ())) ();
    //return (::isInteger (l) && (fabs (l - (*(Ub ())) ()) < COUENNE_EPS));
  }

  /// Set this variable as integer
  virtual inline void setInteger (bool value)
  {integer_ = value ? Integer : Continuous;}

  /// Tell this variable appears once more
  inline void increaseMult () {++multiplicity_;}

  /// Tell this variable appears once less (standardized within
  /// exprSum, for instance)
  inline void decreaseMult () {--multiplicity_;}

  /// Disable this auxiliary variable
  inline void zeroMult () {multiplicity_ = 0;}

  /// How many times this variable appears 
  inline int Multiplicity () {return multiplicity_;}

  /// link this variable to a domain
  inline void linkDomain (Domain *d) {
    domain_ = d; 
    if (lb_) lb_ -> linkDomain (d);
    if (ub_) ub_ -> linkDomain (d);
  }

  /// return proper object to handle expression associated with this
  /// variable (NULL if this is not an auxiliary)
  virtual CouenneObject *properObject (CouenneProblem *p, 
				       Bonmin::BabSetupBase *base, 
				       JnlstPtr jnlst_);
};


/** Structure for comparing expressions
 *
 *  Used in compare() method for same-class expressions
 */

struct compExpr {
  inline bool operator () (exprAux* e0, exprAux* e1) const
  {return (e0 -> Image () -> compare (*(e1 -> Image ())) < 0);}
};


/// allow to draw function within intervals and cuts introduced
void draw_cuts (OsiCuts &, const CouenneCutGenerator *, 
		int, expression *, expression *);

#endif
