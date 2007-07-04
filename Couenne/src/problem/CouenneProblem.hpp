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
#include <expression.hpp>
#include <exprAux.hpp>
#include <CouenneProblemElem.hpp>

struct ASL;
struct expr;


/// structure for comparing expressions

struct compExpr {
  inline bool operator () (exprAux* e0, exprAux* e1) const
  {return (e0 -> Image () -> compare (*(e1 -> Image ())) < 0);}
};


/// The CouenneProblem class

class CouenneProblem {

 protected:

  /// data of the original problem
  std::vector <Objective         *> objectives_;
  std::vector <CouenneConstraint *> constraints_;

  /// data of the linearized problem
  std::vector <exprAux          *> auxiliaries_;
  std::vector <exprVar          *> variables_;

  /// current value and bounds for original- and auxiliary variables
  CouNumber *x_;
  CouNumber *lb_;
  CouNumber *ub_;

  /// expression map for comparison in standardization
  std::set <exprAux *, compExpr> *auxSet_; // int to count occurrences

  /// number of elements in the x_, lb_, ub_ arrays
  int curnvars_;

  /// number of discrete variables
  int nIntVars_;

  /// best solution known (to be loaded from file)
  CouNumber *optimum_;

  /// best known objective function
  CouNumber bestObj_;

  /// indices of variables appearing in products (used for SDP cuts)
  int *quadIndex_;

 public:

  /// constructors, destructor
  CouenneProblem  (const ASL * = NULL);
  CouenneProblem  (const CouenneProblem &);
  ~CouenneProblem ();

  /// clone method (for use within CouenneCutGenerator::clone)
  CouenneProblem *clone () const;

  /// update value of variables and bounds
  void update (CouNumber *, CouNumber *, CouNumber *, int = -1);

  /// get size of vectors
  int nObjs   () const {return objectives_.   size ();}
  int nNLCons () const {return constraints_.  size ();}

  int nAuxs    () const {return auxiliaries_. size ();}
  int nVars    () const {return variables_.   size ();}
  int nIntVars () const {return nIntVars_;}

  /// get elements from vectors
  CouenneConstraint *NLCon (int i) const {return constraints_  [i];}
  Objective         *Obj   (int i) const {return objectives_   [i];}

  ///  LinearConstraint  *Con   (int i) const {return linearconstraints_ [i];}
  exprAux     *Aux   (int i) const {return auxiliaries_  [i];}
  exprVar     *Var   (int i) const {return variables_    [i];}

  /// get and set current variable and bounds
  CouNumber   &X     (int i) {return x_   [i];}
  CouNumber   &Lb    (int i) {return lb_  [i];}
  CouNumber   &Ub    (int i) {return ub_  [i];}

  /// get and set current variable and bounds
  CouNumber  *X     () const {return x_;}
  CouNumber  *Lb    () const {return lb_;}
  CouNumber  *Ub    () const {return ub_;}

  /// get optimal solution and objective value
  CouNumber  *bestSol () const {return optimum_;}
  CouNumber   bestObj () const {return bestObj_;}

  /// add (non linear) objective function
  void addObjective     (expression *, const std::string &);

  /// add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression *);
  void addGEConstraint  (expression *, expression *);
  void addLEConstraint  (expression *, expression *);
  void addRNGConstraint (expression *, expression *, expression *);

  /// add variable
  expression *addVariable (bool);

  /// add auxiliary variable and associate it with expression given as
  /// argument (used in standardization)
  exprAux *addAuxiliary (expression *);

  /// break problem's nonlinear constraints in simple expressions to be
  /// convexified later
  void standardize ();

  /// output objective, linear and nonlinear constraints, and auxiliary
  /// variables
  void print (std::ostream &);

  /// read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL *);

  /// read problem from .nl file using the Ampl Solver Library (ASL)
  expression *nl2e (struct expr *);

  /// store nonlinear problem into a .mod file (with lots of defined
  /// variables)
  ///
  /// @aux controls the use of auxiliaries. If true, a problem is
  /// written with auxiliary variables written with their associated
  /// expression, i.e. w_i = h_i(x,y,w) and bounds l_i <= w_i <= u_i,
  /// while if false these constraints are written in the form l_i <=
  /// h_i (x,y) <= u_i
  /// 
  void writeMod (const std::string &, /// Name of the .mod file to be written
		 bool);               /// true: with aux, false: without.

  /// initialize auxiliary variables and their bounds from original
  /// variables
  void initAuxs (CouNumber *, CouNumber *, CouNumber *);

  /// get auxiliary variables from original variables
  void getAuxs (CouNumber *);

  /// bound tightening
  int tightenBounds (t_chg_bounds *) const;

  /// search for implied bounds 
  int impliedBounds (t_chg_bounds *) const;

  /// look for quadratic terms to be used with SDP cuts
  void fillQuadIndices ();
};


/// read best known solution from file given in argument
bool readOptimum (const std::string &, CouNumber *&, CouNumber &, CouenneProblem *);

#endif
