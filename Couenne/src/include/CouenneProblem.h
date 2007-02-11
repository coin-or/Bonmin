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
struct expr2;

// The CouenneProblem class

class CouenneProblem {

 protected:

  // data of the original problem
  std::vector <Objective         *> objectives_;
  std::vector <CouenneConstraint *> constraints_;

  // data of the linearized problem
  std::vector <exprAux          *> auxiliaries_;
  std::vector <exprVar          *> variables_;

  // current value and bounds for original- and auxiliary variables
  CouNumber *x_;
  CouNumber *lb_;
  CouNumber *ub_;

  // expression map for comparison in standardization
  std::map <std::string, exprAux *> *auxMap_;

 public:

  // constructors, destructor
  CouenneProblem  () {x_ = lb_ = ub_ = NULL; auxMap_ = NULL;}
  CouenneProblem  (const CouenneProblem &);
  ~CouenneProblem ();

  // clone method (for use within CouenneCutGenerator::clone)
  CouenneProblem *clone () const;

  // update value of variables and bounds
  void update (CouNumber *, CouNumber *, CouNumber *);

  // get size of vectors
  int nObjs   () const {return objectives_.        size ();}
  int nNLCons () const {return constraints_.       size ();}

  //  int nCons   () const {return linearconstraints_. size ();}
  int nAuxs   () const {return auxiliaries_.       size ();}
  int nVars   () const {return variables_.         size ();}

  // get elements from vectors
  CouenneConstraint *NLCon (int i) const {return constraints_  [i];}
  Objective         *Obj   (int i) const {return objectives_   [i];}

  //  LinearConstraint  *Con   (int i) const {return linearconstraints_ [i];}
  exprAux     *Aux   (int i) const {return auxiliaries_  [i];}
  exprVar     *Var   (int i) const {return variables_    [i];}

  // get and set current variable and bounds
  CouNumber   &X     (int i) {return x_   [i];}
  CouNumber   &Lb    (int i) {return lb_  [i];}
  CouNumber   &Ub    (int i) {return ub_  [i];}

  // get and set current variable and bounds
  CouNumber  *X     () const {return x_;}
  CouNumber  *Lb    () const {return lb_;}
  CouNumber  *Ub    () const {return ub_;}

  // add (non linear) objective function
  void addObjective     (expression *, const std::string &);

  // add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression *);
  void addGEConstraint  (expression *, expression *);
  void addLEConstraint  (expression *, expression *);
  void addRNGConstraint (expression *, expression *, expression *);

  // add variable
  expression *addVariable (bool);

  // add auxiliary variable and associate it with expression given as
  // argument (used in standardization)
  exprAux *addAuxiliary (expression *);

  // break problem's nonlinear constraints in simple expressions to be
  // convexified later
  void standardize ();

  // output objective, linear and nonlinear constraints, and auxiliary
  // variables
  void print (std::ostream &);

  // read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL_pfgh *);

  // read problem from .nl file using the Ampl Solver Library (ASL)
  expression *nl2e (struct expr2 *);

  // store nonlinear problem into a .mod file (with lots of defined
  // variables)
  void writeMod (char *);
};


// utility to allocate space for coeff, indices, rhs and sign based on
// the data in n and nterms

void allocateCon (int n, int *nterms,                     // input data
 		  expression ***& coeff, int **& indices, // allocated data
		  expression **& rhs, enum con_sign *& sign);
#endif
