/*
 * Name:    expression.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class expression
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRESSION_HPP
#define COUENNE_EXPRESSION_HPP

#include <iostream>
#include <set>
#include <vector>

#include "CouennePrecisions.hpp"
#include "CouenneTypes.hpp"

class OsiBranchingInformation;
class OsiSolverInterface;
class OsiCuts;

class CouenneProblem;
class CouenneCutGenerator;
class CouenneObject;

class exprAux;
class exprCopy;
class exprVar;

class DepNode;
class DepGraph;

class Domain;

struct compNode;

/// Expression base class
///
/// An empty expression class with no type or operator() from which
/// all other expression classes (for constants, variables, and
/// operators) are derived.

class expression {

 public:

  /// Constructor
  expression () {}

  /// Copy constructor. Pass pointer to variable vector when
  /// generating new problem, whose set of variables is equivalent but
  /// may be changed or whose value is independent.
  expression (const expression &e, Domain *d = NULL) {}

  /// Destructor
  virtual ~expression () {}

  /// Cloning method
  virtual expression *clone (Domain *d = NULL) const 
  {return NULL;}

  /// Return index of variable (only valid for exprVar and exprAux)
  virtual inline int Index () const
  {return -1;}

  /// return number of arguments (when applicable, that is, with N-ary functions)
  virtual inline int nArgs () const
  {return 0;}

  /// return arglist (when applicable, that is, with N-ary functions)
  virtual inline expression **ArgList () const
  {return NULL;}

  /// set arglist (used in deleting nodes without deleting children)
  virtual inline void ArgList (expression **al)
  {}

  /// return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
  {return NULL;}

  /// return pointer to argument (when applicable, i.e., with univariate functions)
  virtual inline expression **ArgPtr ()
  {return NULL;}

  /// node type
  virtual inline enum nodeType Type () const
  {return EMPTY;}

  /// return pointer to corresponding expression (for auxiliary variables only)
  virtual inline expression *Image () const
  {return NULL;}

  /// set expression associated with this auxiliary variable (for
  /// compatibility with exprAux)
  void Image (expression *image) {}

  /// value (empty)
  virtual inline CouNumber Value () const 
  {return 0.;}

  /// If this is an exprClone of a exprClone of an expr???, point to
  /// the original expr??? instead of an exprClone -- improve computing
  /// efficiency. Only overloaded for exprClones/exprCopy, of course.
  virtual inline const expression *Original () const 
    {return this;}

  /// print expression to iostream
  virtual void print (std::ostream &s  = std::cout,   /// output stream
		      bool             = false) const /// descend into auxiliaries' image?
    {s << '?';}

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
  {return 0;}

  /// differentiation
  virtual inline expression *differentiate (int) 
  {return NULL;}

  /// dependence on variable set: return cardinality of subset of the
  /// set of indices in first argument which occur in expression. 
  virtual int dependsOn (int *ind, int n, 
			 enum dig_type type = STOP_AT_AUX);

  /// version with one index only
  inline int dependsOn (int singleton, 
			 enum dig_type type = STOP_AT_AUX)
  {return dependsOn (&singleton, 1, type);}

  /// fill std::set with indices of variables on which this expression
  /// depends. Also deal with expressions that have no variable
  /// pointers (exprGroup, exprQuad)
  virtual inline int DepList (std::set <int> &deplist, 
			      enum dig_type   type = ORIG_ONLY)
  {return 0;}

  /// simplify expression (useful for derivatives)
  virtual inline expression *simplify () 
  {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
  {return NONLINEAR;}

  /// is this expression defined as an integer?
  virtual inline bool isDefinedInteger ()
  {return isInteger ();}

  /// is this expression integer?
  virtual inline bool isInteger ()
  {return false;}

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Get lower and upper bound of an expression (if any) -- real values
  virtual void getBounds (CouNumber &, CouNumber &);

  /// Create standard form of this expression, by:
  ///
  /// - creating auxiliary w variables and corresponding expressions
  /// - returning linear counterpart as new constraint (to replace 
  ///   current one)
  ///
  /// For the base exprOp class we only do the first part (for argument
  /// list components only), and the calling class (Sum, Sub, Mul, Pow,
  /// and the like) will do the part for its own object
  ///
  /// addAux is true if a new auxiliary variable should be added
  /// associated with the standardized expression
  virtual inline exprAux *standardize (CouenneProblem *p, bool addAux = true) 
  {return NULL;}

  /// generate convexification cut for constraint w = this
  virtual void generateCuts (expression *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg = NULL, int wind = -1, 
			     CouNumber lb = -COUENNE_INFINITY, 
			     CouNumber ub =  COUENNE_INFINITY) {}

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ()
  {return NULL;}

  /// return integer for comparing expressions (used to recognize
  /// common expression)
  virtual enum expr_type code () 
  {return COU_EXPRESSION;}

  /// either CONVEX, CONCAVE, AFFINE, or NONCONVEX
  virtual enum convexity convexity () 
  {return NONCONVEX;}

  /// compare expressions
  virtual int compare (expression &);

  /// compare copies of expressions
  virtual int compare (exprCopy   &);

  /// used in rank-based branching variable choice: original variables
  /// have rank 1; auxiliary w=f(x) has rank r(w) = r(x)+1; finally,
  /// auxiliary w=f(x1,x2...,xk) has rank r(w) = 1+max{r(xi):i=1..k}.
  virtual int rank ()
  {return -1;} // return null rank

  /// does a backward implied bound processing on every expression,
  /// including exprSums although already done by Clp (useful when
  /// repeated within Couenne). Parameters are the index of the
  /// (auxiliary) variable in question and the current lower/upper
  /// bound. The method returns true if there has been a change on any
  /// bound on the variables on which the expression depends.
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *)
  {return false;}

  /// multiplicity of a variable
  virtual inline int Multiplicity ()
  {return 1;}

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments. Return estimated improvement in
  /// objective function
  virtual CouNumber selectBranch (const CouenneObject *obj, 
				  const OsiBranchingInformation *info,
				  expression * &var, 
				  double * &brpts, 
				  int &way)
  {var = NULL; return 0.;}

  /// replace expression with another
  virtual void replace (exprVar *, exprVar *) {}

  /// update dependence set with index of variables on which this
  /// expression depends
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *) {}

  /// return pointer to variable domain
  virtual Domain *domain () 
  {return NULL;}

  /// empty function to update domain pointer
  virtual void linkDomain (Domain *d) {}

  /// empty function to redirect variables to proper variable vector
  virtual void realign (const CouenneProblem *p) {}

  /// indicating if function is monotonotically increasing
  virtual bool isBijective() const {return false};

  /// compute the inverse function
  virtual CouNumber inverse(expression *vardep) const {
    return -COUENNE_INFINITY;
  }

  /// closest feasible points in function in both directions
  virtual void closestFeasible (expression *varind, expression *vardep,
				CouNumber& left, CouNumber& right)
  {
    assert(isBijective());
    CouNumber inv = inverse(vardep);
    CouNumber curr = (*varind) ();
    if (curr > inv) {
      left  = inv;
      right = curr;
    }
    else {
      left  = curr;
      right = inv;
    }
  }

};


/// updates maximum violation. Used with all impliedBound. Returns true
/// if a bound has been modified, false otherwise

inline bool updateBound (int sign, CouNumber *dst, CouNumber src) {

  // meaning: 
  //
  // if (*dst > src) && (sign > 0) --> dst down to src
  // if (*dst < src) && (sign < 0) --> dst up   to src
  //
  // that is, sign > 0 means we are tightening an UPPER bound
  //          sign < 0                            LOWER

  if (((sign > 0) ? (*dst - src) : (src - *dst)) > COUENNE_EPS) {
    //printf ("%.12g --> %.12g\n", *dst, src);
    *dst = src; // tighten
    return true;
  }

  return false;
}

/// independent comparison
inline int compareExpr (const void *e0, const void *e1) {
  return ((*(expression **) e0) -> compare (**(expression **)e1));
}

/// is this number integer?
inline bool isInteger (CouNumber x)
{return (fabs (COUENNE_round (x) - x) < COUENNE_EPS_INT);}

#endif
