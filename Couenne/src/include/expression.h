/*
 * Name:    expression.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class expression
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRESSION_H
#define COUENNE_EXPRESSION_H

#define STACK_SIZE 10000

#include <iostream>
#include <CouennePrecisions.h>
#include <CouenneTypes.h>

class CouenneProblem;
class CouenneCutGenerator;
class exprAux;
class OsiSolverInterface;
class OsiCuts;


// expression base class

class expression {

 protected:

  // Static members to be used "globally" by an expression upon
  // evaluation with a given value of the variables' or the bounds'
  // vectors. 
  //
  // stack is used in evaluation as a LIFO structure where all
  // arguments of an expression (which are expressions themselves) are
  // stored (with a PUSH-like operation) after being evaluated, for
  // the current node to process them with a POP operation. PUSH and
  // POP operations are the *++sp and *sp-- instructions,
  // respectively, on the Stack Pointer variable sp.
  // 
  // STACK_SIZE should be enough for expressions where at each node of
  // the evaluation tree, DEPTH + #ARGUMENTS is at most STACK_SIZE,
  // where DEPTH is the depth of the evaluation node and #ARGUMENTS is
  // the number of arguments of the function in the node.

  static CouNumber stack [STACK_SIZE];
  static CouNumber *sp;

  // these "global" variables, static members of expression, contain
  // the current value of the variables' and bounds' vectors. The
  // former vector is used to compute the value of the expression, for
  // instance when using a non-linear solver that requires evaluation
  // of all expressions in the problem. The latter is used when
  // updating the problem after a change in the variables' bounds.
  //
  // CAUTION: every time an expression (or a set) is evaluated, the
  // user can (SHOULD) tell the program to synchronize with her/his
  // own vectors by calling the method expression::update (myvars,
  // mylbounds, myubounds) below.

  static CouNumber *variables_;
  static CouNumber *lbounds_;
  static CouNumber *ubounds_;

  // current value of the expression, used when accessing a copy of
  // the expression created for a node that is evaluated after the
  // original (saves some time).

  CouNumber currValue_;

 public:

  // update the value of "external" vectors (if the addresses have not
  // changed, it is not needed)
  static inline void update (CouNumber *variables = NULL, 
			     CouNumber *lbounds   = NULL, 
			     CouNumber *ubounds   = NULL) {

    if (variables) variables_ = variables;
    if (lbounds)   lbounds_   = lbounds;
    if (ubounds)   ubounds_   = ubounds;
  }

  // Constructor, destructor
  expression () {}
  expression (const expression &e) {}
  virtual ~expression () {}

  // cloning method
  virtual expression *clone () {return NULL;}

  // return index of variable (only valid for exprVar and exprAux)
  virtual inline int Index () const
  {return -1;}

  // return arglist (when applicable, that is, with N-ary functions)
  virtual inline expression **ArgList () const
    {return NULL;}

  // return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
    {return NULL;}

  // node type
  virtual inline enum nodeType Type ()
    {return EMPTY;}

  // value (empty)
  virtual inline CouNumber Value () const 
    {return currValue_;}

  // If this is an exprClone of a exprClone of an expr???, point to
  // the original expr??? instead of an exprClone -- improve computing
  // efficiency. Only overloaded for exprClones/exprCopy, of course.
  virtual inline expression *Original () {return this;}

  // String equivalent (for comparisons)
  virtual std::string name() {return "";}

  // I/O
  virtual void print (std::ostream &) {}

  // null function for evaluating the expression
  virtual inline CouNumber operator () () 
    {return 0;}

  // differentiation
  virtual inline expression *differentiate (int) 
    {return NULL;}

  // dependence on variable set
  virtual inline bool dependsOn (int *, int) 
    {return false;}

  // simplify expression (useful for derivatives)
  virtual inline expression *simplify () 
    {return NULL;}

  // get a measure of "how linear" the expression is:
  //
  // 0: a constant
  // 1: linear
  // 2: quadratic
  // 3: nonlinear non-quadratic
  virtual inline int Linearity ()
    {return NONLINEAR;}

  // Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //
  // For convex functions $w = f(x)$, add linear inequality
  //
  // $w - f'(x_k) x   \ge   f(x_k) - f' (x_k) x_k$, or
  // $lhs$            \ge   $rhs$
  //
  // where $w$ is the aux variable, $x$ is the unique variable of the
  // term, $f'$ is the derivative with respect to $x$, and all $x_k$ are
  // all the nPWApprox() sample points (at least two: lower- and upper
  // bound of $x$

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  /*  virtual inline int lowerLinearHull (exprAux *, int *&, expression ***&, 
			       int **&, expression **&, enum con_sign *&)
    {return 0;}

  // similarly, construct linear over-estimator for expression within
  // problem *p (p is used to add convexification constraints). It is
  // also used when this function appears with a minus sign in the
  // expression
  virtual inline int upperLinearHull (exprAux *, int *&, expression ***&, 
			       int **&, expression **&, enum con_sign *&)
    {return 0;}
  */
  // Create standard formulation of this expression, by:
  //
  // - creating auxiliary w variables and corresponding expressions
  // - returning linear counterpart as new constraint (to replace 
  //   current one)
  //
  // For the base exprOp class we only do the first part (for argument
  // list components only), and the calling class (Sum, Sub, Mul, Pow,
  // and the like) will do the part for its own object
  virtual inline exprAux *standardize (CouenneProblem *) 
    {return NULL;}

  // generate convexification cut for constraint w = this
  virtual void generateCuts (exprAux *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *) {}
};

#endif
