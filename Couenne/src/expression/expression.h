/*
 * Name:    expression.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class expression
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRESSION_H
#define COUENNE_EXPRESSION_H

#define STACK_SIZE 100000

#include <iostream>

#include <CouennePrecisions.h>
#include <CouenneTypes.h>

class OsiBranchingInformation;
class OsiBranchingObject;
class CouenneProblem;
class CouenneCutGenerator;
class exprAux;
class OsiSolverInterface;
class OsiCuts;
class exprUnary;
class exprOp;
class exprCopy;

/// expression base class

class expression {

 protected:

  /// Static members to be used "globally" by an expression upon
  /// evaluation with a given value of the variables' or the bounds'
  /// vectors. 
  ///
  /// stack is used in evaluation as a LIFO structure where all
  /// arguments of an expression (which are expressions themselves) are
  /// stored (with a PUSH-like operation) after being evaluated, for
  /// the current node to process them with a POP operation. PUSH and
  /// POP operations are the *++sp and *sp-- instructions,
  /// respectively, on the Stack Pointer variable sp.
  /// 
  /// STACK_SIZE should be enough for expressions where at each node of
  /// the evaluation tree, DEPTH + #ARGUMENTS is at most STACK_SIZE,
  /// where DEPTH is the depth of the evaluation node and #ARGUMENTS is
  /// the number of arguments of the function in the node.
  ///
  /// UPDATE (04/12/07): no longer used, keep for possible new
  /// operators

  static CouNumber stack [STACK_SIZE];
  static CouNumber *sp;

  /// these "global" variables, static members of expression, contain
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
  static CouNumber *lbounds_;
  static CouNumber *ubounds_;

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

  /// return current values of variables and bounds
  static CouNumber Lbound   (int i) {return lbounds_   [i];}
  static CouNumber Ubound   (int i) {return ubounds_   [i];}
  static CouNumber Variable (int i) {return variables_ [i];}

  /// return arrays
  static CouNumber *Lbounds   () {return lbounds_;}
  static CouNumber *Ubounds   () {return ubounds_;}
  static CouNumber *Variables () {return variables_;}

  /// Constructor, destructor
  expression () {}
  expression (const expression &e) {}
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

  /// I/O
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

  /// dependence on variable set
  virtual inline bool dependsOn (int *, int) 
    {return false;}

  /// simplify expression (useful for derivatives)
  virtual inline expression *simplify () 
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return NONLINEAR;}

  /// is this expression integer?
  virtual bool isInteger ()
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
  virtual inline exprAux *standardize (CouenneProblem *) 
    {return NULL;}

  /// generate convexification cut for constraint w = this
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL) {}

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
  /// bound. The method returns the best bound improvement obtained on
  /// all variables of the expression.
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

  /// distance covered by current point if branching rule applied to this expression
  /*virtual double BranchGain (expression *, const OsiBranchingInformation *)
    {return 0.;}

  /// branching object best suited for this expression
  virtual OsiBranchingObject *BranchObject (expression *, const OsiBranchingInformation *)
  {return NULL;}*/
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
    //printf ("%.9f --> %.9f\n", *dst, src);
    *dst = src; // tighten
    return true;
  }

  return false;
}

/// independent comparison
inline int compareExpr (const void *e0, const void *e1) {
  return ((*(expression **) e0) -> compare (**(expression **)e1));
}

/// maximum
inline CouNumber mymin (register CouNumber a, register CouNumber b) 
{return ((a<b) ? a : b);} 

/// minimum
inline CouNumber mymax (register CouNumber a, register CouNumber b) 
{return ((a>b) ? a : b);} 

#endif
