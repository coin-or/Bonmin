/*
 * Name:    CouenneProblem.h
 * Author:  Pietro Belotti
 * Purpose: define the class CouenneProblem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PROBLEM_H
#define COUENNE_PROBLEM_H

#include <vector>

#include <OsiRowCut.hpp>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <CouenneProblemElem.h>

struct ASL_pfgh;

// The CouenneProblem class

class CouenneProblem {

 protected:

  // data of the original problem
  std::vector <Objective         *> objectives_;
  std::vector <CouenneConstraint *> constraints_;

  // data of the linearized problem
  std::vector <LinearConstraint *> linearconstraints_;
  std::vector <exprAux          *> auxiliaries_;
  std::vector <exprVar          *> variables_;

  // current value and bounds for original- and auxiliary variables
  CouNumber *x_;
  CouNumber *lb_;
  CouNumber *ub_;

 public:

  // constructor, destructor
  CouenneProblem  () {}
  ~CouenneProblem ();

  // update value of variables and bounds
  void update (CouNumber *x, CouNumber *l, CouNumber *u) 
  {x_ = x; lb_ = l; ub_ = u;}

  // get size of vectors
  int nObjs   () const {return objectives_.        size ();}
  int nNLCons () const {return constraints_.       size ();}

  int nCons   () const {return linearconstraints_. size ();}
  int nAuxs   () const {return auxiliaries_.       size ();}
  int nVars   () const {return variables_.         size ();}

  // get elements from vectors
  CouenneConstraint *NLCon (int i) const {return constraints_       [i];}
  Objective         *Obj   (int i) const {return objectives_        [i];}

  LinearConstraint  *Con   (int i) const {return linearconstraints_ [i];}
  exprAux           *Aux   (int i) const {return auxiliaries_       [i];}
  exprVar           *Var   (int i) const {return variables_         [i];}

  // get and set current variable and bounds
  CouNumber         &X     (int i) {return x_                 [i];}
  CouNumber         &Lb    (int i) {return lb_                [i];}
  CouNumber         &Ub    (int i) {return ub_                [i];}

  // get and set current variable and bounds
  CouNumber         *X     () {return x_;}
  CouNumber         *Lb    () {return lb_;}
  CouNumber         *Ub    () {return ub_;}

  // add (non linear) objective function
  void addObjective     (expression *, const std::string &);

  // add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression *);
  void addGEConstraint  (expression *, expression *);
  void addLEConstraint  (expression *, expression *);
  void addRNGConstraint (expression *, expression *, expression *);

  // add linear constraints defined as - w + a * x >=< b, where a is
  // an expression (function of variable bounds) pointed to by coeff,
  // where b is the expression rhs
  void inline addLinearConstraint (exprAux *w, exprVar *x, expression *coeff, 
				   expression *rhs, enum con_sign sign)
    {linearconstraints_ . push_back (new LinearConstraint  (w, x, coeff, rhs, sign));}

  // add linear constraint defined by array of n coefficients and
  // array of n indices, with lower- and upper bounds and sign
  void inline addLinearConstraint (int n, expression **coeffs, int *indices,
				   expression *rhs, enum con_sign sign)
    {linearconstraints_ . push_back (new LinearConstraint (n, coeffs, indices, rhs, sign));}

  // add variable
  expression *addVariable ();

  // add auxiliary variable and associate it with expression given as
  // argument (used in standardization)
  exprAux *addAuxiliary (expression *);

  // break problem's nonlinear constraints in simple expressions to be
  // convexified later
  void standardize ();

  // convexify standardized problem
  //  void convexify ();

  // output objective, linear and nonlinear constraints, and auxiliary
  // variables
  void print (std::ostream &);

  // read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL_pfgh *);
};


// utility to allocate space for coeff, indices, rhs and sign based on
// the data in n and nterms

void allocateCon (int n, int *nterms,                     // input data
 		  expression ***& coeff, int **& indices, // allocated data
		  expression **& rhs, enum con_sign *& sign);
#endif
