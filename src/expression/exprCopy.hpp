/*
 * Name:    exprCopy.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprCopy
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOPY_HPP
#define COUENNE_EXPRCOPY_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"

class CouenneObject;

// expression copy (points to VALUE of another expression) 

class exprCopy: public expression {

 protected:

  /// the expression this object is a (reference) copy of
  expression *copy_;

  /// saved value to be used by exprStore expressions
  CouNumber value_;

 public:

  /// node type
  inline enum nodeType Type () const
    {return copy_ -> Type ();}

  /// Constructor
  exprCopy  (expression *copy):
    copy_  (copy),
    value_ (0.) {}

  /// Copy constructor
  exprCopy (const exprCopy &e, Domain *d = NULL);

  /// Destructor -- CAUTION: this is the only destructive destructor,
  /// exprClone and exprStore do not destroy anything
  virtual ~exprCopy () {
    if (copy_)
      delete copy_;
  }

  /// Cloning method
  virtual inline expression *clone (Domain *d = NULL) const
  {return new exprCopy (*this, d);}

  /// If this is an exprClone of a exprClone of an expr???, point to
  /// the original expr??? instead of an exprClone -- improves computing
  /// efficiency
  inline const expression *Original () const
  {return copy_ -> Original ();}

  /// return pointer to corresponding expression (for auxiliary variables only)
  virtual inline expression *Image () const
  {return copy_ -> Image ();}

  /// Get variable index in problem
  inline int Index () const
  {return copy_ -> Index ();}

  /// Return number of arguments (when applicable, that is, with N-ary functions)
  virtual inline int nArgs () const
  {return copy_ -> nArgs ();}

  /// return arglist (when applicable, that is, with N-ary functions)
  virtual inline expression **ArgList () const
  {return copy_ -> ArgList ();}

  /// set arglist (used in deleting nodes without deleting children)
  virtual inline void ArgList (expression **al)
  {copy_ -> ArgList (al);}

  /// return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
  {return copy_ -> Argument ();}

  /// return pointer to argument (when applicable, i.e., with univariate functions)
  virtual inline expression **ArgPtr ()
  {return copy_ -> ArgPtr ();}

  /// I/O
  virtual void print (std::ostream &out = std::cout, 
		      bool descend      = false) const
  {copy_ -> Original () -> print (out, descend);}
  //{out << "["; copy_ -> print (out, descend); out << "]"; } // Must go

  /// value
  virtual inline CouNumber Value () const 
  {return value_;}

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
  {return (value_ = (*copy_) ());}
    //    {return (*copy_) ();}
    //    {return (copy_ -> Value ());}

  /// return l-2 norm of gradient at given point
  virtual inline CouNumber gradientNorm (const double *x)
  {return copy_ -> gradientNorm (x);}

  /// differentiation
  inline expression *differentiate (int index) 
  {return copy_ -> differentiate (index);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  inline int DepList (std::set <int> &deplist, 
		      enum dig_type   type = ORIG_ONLY)
  {return copy_ -> DepList (deplist, type);}

  /// simplify expression (useful for derivatives)
  inline expression *simplify () 
  {return copy_ -> simplify ();}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
  {return copy_ -> Linearity ();}

  virtual inline bool isInteger ()
  {return copy_ -> isInteger ();}

  /// Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) 
  {copy_ -> getBounds (lower, upper);}

  /// Create standard formulation of this expression
  inline exprAux *standardize (CouenneProblem *p, bool addAux = true)
  {return copy_ -> standardize (p, addAux);}

  /// generate convexification cut for constraint w = this
  inline void generateCuts (expression *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg = NULL, int wind= -1, 
			    CouNumber lb = -COUENNE_INFINITY, 
			    CouNumber ub =  COUENNE_INFINITY)

  {copy_ -> generateCuts (w, si, cs, cg, chg, wind, lb, ub);}

  /// code for comparisons
  inline enum expr_type code () 
  {return copy_ -> code ();}

  /// either CONVEX, CONCAVE, AFFINE, or NONCONVEX 
  virtual inline enum convexity convexity ()
  {return copy_ -> convexity ();}

  /// compare this with other expression
  int compare (expression &e) 
  {return copy_ -> compare (e);}

  /// used in rank-based branching variable choice
  inline int rank ()
  {return copy_ -> rank ();} 

  /// implied bound processing
  inline bool impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg)
  {return copy_ -> impliedBound (wind, l, u, chg);}

  /// multiplicity of a variable: how many times this variable occurs
  /// in expressions throughout the problem
  virtual inline int Multiplicity ()
  {return copy_ -> Multiplicity ();}

  /// Set up branching object by evaluating many branching points for each expression's arguments.
  /// Return estimated improvement in objective function 
  virtual inline CouNumber selectBranch (const CouenneObject *obj, 
				  const OsiBranchingInformation *info,
				  expression * &var, 
				  double * &brpts, 
				  double * &brDist, // distance of current LP
						    // point to new convexifications
				  int &way)

  {return copy_ -> selectBranch (obj, info, var, brpts, brDist, way);}

  /// replace occurrence of a variable with another variable
  virtual void replace (exprVar *, exprVar *);

  /// fill in dependence structure
  inline void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g)
  {copy_ -> fillDepSet (dep, g);}

  /// redirect variables to proper variable vector
  virtual void realign (const CouenneProblem *p)
  {copy_ -> realign (p);}

  /// indicating if function is monotonically increasing
  virtual bool isBijective() const 
  {return copy_ -> isBijective ();}

  /// compute the inverse function
  virtual CouNumber inverse (expression *vardep) const
  {return copy_ -> inverse (vardep);}

  /// closest feasible points in function in both directions
  virtual void closestFeasible (expression *varind, expression *vardep,
				CouNumber& left, CouNumber& right) const
  {copy_ -> closestFeasible (varind, vardep, left, right);}
};

#endif
