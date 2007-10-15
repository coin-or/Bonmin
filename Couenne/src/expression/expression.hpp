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

#define STACK_SIZE 10000

#include <iostream>
#include <set>
#include <cmath>
#include <algorithm>
#include <iterator>

#include <CouennePrecisions.hpp>
#include <CouenneTypes.hpp>

class OsiBranchingInformation;
class OsiBranchingObject;
class OsiSolverInterface;
class OsiCuts;

class CouenneProblem;
class CouenneCutGenerator;

class exprAux;
class exprUnary;
class exprOp;
class exprCopy;
class exprVar;

class DepNode;
class DepGraph;

struct compNode;

/// Expression base class
///
/// An empty expression class with no type or operator() from which
/// all other expression classes (for constants, variables, and
/// operators) are derived.

class expression {

 protected:

  /// Static members to be used "globally" by an expression upon
  /// evaluation with a given value of the variables' or the bounds'
  /// vectors.
  ///
  /// stack is used in evaluation as a LIFO structure where all
  /// arguments of an expression (which are expressions themselves)
  /// are stored (with a PUSH operation) after being evaluated, for
  /// the current node to process them with a POP operation. PUSH and
  /// POP operations are the *++sp and *sp-- instructions,
  /// respectively, on the Stack Pointer variable sp.
  /// 
  /// STACK_SIZE should be enough for expressions where at each node of
  /// the evaluation tree, DEPTH + #ARGUMENTS is at most STACK_SIZE,
  /// where DEPTH is the depth of the evaluation node and #ARGUMENTS is
  /// the number of arguments of the function in the node.

  static CouNumber stack [STACK_SIZE];

  /// Stack pointer: a cursor on a LIFO structure to hold the recently
  /// computed argument(s) of an expression

  static CouNumber *sp;

  /// These "global" variables, static members of expression, contain
  /// the current value of the variables' and bounds' vectors. The
  /// former vector is used to compute the value of the expression, for
  /// instance when using a non-linear solver that requires evaluation
  /// of all expressions in the problem. The latter is used when
  /// updating the problem after a change in the variables' bounds.
  ///
  /// CAUTION: every time an expression (or a set) is evaluated, the
  /// user can (SHOULD) tell the program to synchronize with her/his
  /// own vectors by calling the method expression::update (myvars,
  /// mylbounds, myubounds) below.

  static CouNumber *variables_;
  static CouNumber *lbounds_;   ///< vector of lower bounds
  static CouNumber *ubounds_;   ///< vector of upper bounds

  /// current value of the expression, used when accessing a copy of
  /// the expression created for a node that is evaluated after the
  /// original (saves some time).

  CouNumber currValue_;

 public:

  /// update the value of "external" vectors (if the addresses have not
  /// changed, it is not needed)
  static inline void update (CouNumber *variables = NULL, 
			     CouNumber *lbounds   = NULL, 
			     CouNumber *ubounds   = NULL) {

    if (variables) variables_ = variables;
    if (lbounds)   lbounds_   = lbounds;
    if (ubounds)   ubounds_   = ubounds;
  }

  // return current values of variables and bounds
  static CouNumber Lbound   (int i) {return lbounds_   [i];} ///< return \f$l_i\f$
  static CouNumber Ubound   (int i) {return ubounds_   [i];} ///< return \f$u_i\f$
  static CouNumber Variable (int i) {return variables_ [i];} ///< return \f$x_i\f$

  // return whole vectors
  static CouNumber *Lbounds   () {return lbounds_;}   ///< return vector of lower bounds
  static CouNumber *Ubounds   () {return ubounds_;}   ///< return vector of upper bounds
  static CouNumber *Variables () {return variables_;} ///< return vector of variables

  /// Constructor
  expression () {}

  /// Copy constructor
  expression (const expression &e) {}

  /// Destructor
  virtual ~expression () {}

  /// cloning method
  virtual expression *clone () const 
    {return new expression (*this);}

  /// return index of variable (only valid for exprVar and exprAux)
  virtual inline int Index () const
    {return -1;}

  /// return number of arguments (when applicable, that is, with N-ary functions)
  virtual inline int nArgs () const
    {return 0;}

  /// return arglist (when applicable, that is, with N-ary functions)
  virtual inline expression **ArgList () const
    {return NULL;}

  /// return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
    {return NULL;}

  /// return pointer to argument (when applicable, i.e., with univariate functions)
  virtual inline expression **ArgPtr ()
    {return NULL;}

  /// node type
  virtual inline enum nodeType Type ()
    {return EMPTY;}

  /// value (empty)
  virtual inline CouNumber Value () const 
    {return currValue_;}

  /// If this is an exprClone of a exprClone of an expr???, point to
  /// the original expr??? instead of an exprClone -- improve computing
  /// efficiency. Only overloaded for exprClones/exprCopy, of course.
  virtual inline const expression *Original () const 
    {return this;}

  /// print expression to iostream
  virtual void print (std::ostream &s = std::cout,    /// output stream
		      bool = false,                   /// descend into auxiliaries' image?
		      CouenneProblem * = NULL) const  /// problem pointer (in exprGroup)
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
			 CouenneProblem *p = NULL, 
			 enum dig_type   type = STOP_AT_AUX);

  /// fill std::set with indices of variables on which this expression
  /// depends. Also deal with expressions that have no variable
  /// pointers (exprGroup, exprQuad)
  virtual inline int DepList (std::set <int> &deplist, 
			      enum dig_type   type = ORIG_ONLY,
			      CouenneProblem *p    = NULL)
    {return 0;}

  /// simplify expression (useful for derivatives)
  virtual inline expression *simplify () 
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return NONLINEAR;}

  /// is this expression integer?
  virtual inline bool isInteger ()
    {return false;}

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

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
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg,
			     t_chg_bounds *chg = NULL, int wind = -1, 
			     CouNumber lb = -COUENNE_INFINITY, 
			     CouNumber ub =  COUENNE_INFINITY) {}

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ()
    {printf ("Warning: expression::getFixIndex()\n"); return NULL;}

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
  virtual int rank (CouenneProblem *p = NULL)
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
  virtual int Multiplicity () 
    {return 1;}

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments. Return estimated improvement in
  /// objective function
  virtual CouNumber selectBranch (expression *w, 
				  const OsiBranchingInformation *info,
				  int &ind, 
				  double * &brpts, 
				  int &way)
    {ind = -1; return 0.;}

  /// replace expression with another
  virtual void replace (exprVar *, exprVar *) {}

  /// update dependence set with index of variables on which this
  /// expression depends
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *) {}
};


/// updates maximum violation. Used with all impliedBound. Returns true
/// if a bound has been modified, false otherwise

inline bool updateBound (int sign, CouNumber *dst, CouNumber src) {

  register CouNumber delta = src - *dst;

  if (sign > 0) 
    delta = - delta;

  // meaning: 
  //
  // if (*dst > src) && (sign > 0) --> dst down to src
  // if (*dst < src) && (sign < 0) --> dst up   to src
  //
  // that is, sign > 0 means we are tightening an UPPER bound
  //          sign < 0                            LOWER

  if (delta > COUENNE_EPS) {
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

#endif
