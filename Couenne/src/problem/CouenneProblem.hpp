/*
 * Name:    CouenneProblem.hpp
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

class DepGraph;


/** Structure for comparing expressions
 *
 *  Used in compare() method for same-class expressions
 */

struct compExpr {
  inline bool operator () (exprAux* e0, exprAux* e1) const
  {return (e0 -> Image () -> compare (*(e1 -> Image ())) < 0);}
};


/** Class for MINLP problems
 *
 *  It is read from an AMPL .nl file and contains variables, AMPL's
 *  "defined variables" (aka common expressions), objective(s), and
 *  constraints in the form of expression's. Changes throughout the
 *  program occur in standardization.
 */

class CouenneProblem {

 protected:

  std::vector <CouenneObjective  *> objectives_;  ///< objectives
  std::vector <CouenneConstraint *> constraints_; ///< constraints
  std::vector <exprVar           *> variables_;   ///< variables (original, auxiliary, and defined)

  /// AMPL's common expressions (read from AMPL through structures cexps and cexps1)
  std::vector <expression *> commonexprs_; 

  CouNumber *x_;  ///< current value of all variables
  CouNumber *lb_; ///< current lower bound on all variables
  CouNumber *ub_; ///< current upper bound on all variables

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

  /// variables that have commuted to auxiliary
  bool *commuted_;

  /// numbering of variables. No variable xi with associated pi(i)
  /// greater than pi(j) should be evaluated before variable xj
  //int *numbering_;

  /// number of "defined variables" (aka "common expressions")
  int ndefined_;

  /// dependence (acyclic) graph: shows dependence of all auxiliary
  /// variables on one another and on original variables. Used to
  /// create a numbering of all variables for evaluation and bound
  /// tightening (reverse order for implied bounds)
  DepGraph *graph_;

  /// number of original variables
  int nOrig_;

 public:

  CouenneProblem  (const ASL * = NULL);     ///< Constructor
  CouenneProblem  (const CouenneProblem &); ///< Copy constructor
  ~CouenneProblem ();                       ///< Destructor

  /// clone method (for use within CouenneCutGenerator::clone)
  CouenneProblem *clone () const;

  /// update value of variables and bounds
  void update (CouNumber *, CouNumber *, CouNumber *, int = -1);

  int nObjs   () const {return objectives_.   size ();} ///< get number of objectives
  int nNLCons () const {return constraints_.  size ();} ///< get number of constraints


  int nOrig    () const {return nOrig_;}              ///< number of original (independent) variables
  int nIntVars () const {return nIntVars_;}           ///< number of original integer variables
  int nVars    () const {return variables_. size ();} ///< total number of variables

  // get elements from vectors
  CouenneConstraint *NLCon (int i) const {return constraints_  [i];} ///< pointer to i-th constraint
  CouenneObjective  *Obj   (int i) const {return objectives_   [i];} ///< pointer to i-th objective

  /// return i-th variable
  exprVar *Var   (int i) const {return variables_ [i];}

  /// return vector of variables (symbolic representation
  std::vector <exprVar *> &Variables   () {return variables_;}

  /// return pointer to set for comparisons
  std::set <exprAux *, compExpr> *AuxSet () {return auxSet_;}

  /// return pointer to dependence graph
  DepGraph *getDepGraph () {return graph_;}

  // get and set current variable and bounds
  CouNumber   &X     (int i) {return x_   [i];} ///< return current \f$x_i\f$
  CouNumber   &Lb    (int i) {return lb_  [i];} ///< return current lower bound on \f$x_i\f$
  CouNumber   &Ub    (int i) {return ub_  [i];} ///< return current upper bound on\f$x_i\f$

  // get and set current variable and bounds
  CouNumber  *X     () const {return x_;}  ///< return vector of variables
  CouNumber  *Lb    () const {return lb_;} ///< return vector of lower bounds
  CouNumber  *Ub    () const {return ub_;} ///< return vector of upper bounds

  // get optimal solution and objective value
  CouNumber  *bestSol () const {return optimum_;} ///< best known solution (read from file)
  CouNumber   bestObj () const {return bestObj_;} ///< objective of best known solution

  /// get vector of commuted variables
  bool *&Commuted () {return commuted_;}

  /// add (non linear) objective function
  void addObjective     (expression *, const std::string &);

  // add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression *); ///< add equality (\f$ h(x) = b\f$) constraint
  void addGEConstraint  (expression *, expression *); ///< add  (\f$h(x)\ge b\f$) constraint
  void addLEConstraint  (expression *, expression *); ///< add  (\f$h(x)\le b\f$) constraint
  void addRNGConstraint (expression *, expression *, 
			 expression *);               ///< add range (\f$a\le h(x)\le b\f$) constraint

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
  void print (std::ostream & = std::cout);

  /// read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL *);

  /// read problem from .nl file using the Ampl Solver Library (ASL)
  expression *nl2e (struct expr *);

  /// store nonlinear problem into a .mod file (with lots of defined
  /// variables)
  /// 
  /// @param fname Name of the .mod file to be written
  ///
  /// @param aux controls the use of auxiliaries. If true, a problem
  /// is written with auxiliary variables written with their
  /// associated expression, i.e. w_i = h_i(x,y,w) and bounds l_i <=
  /// w_i <= u_i, while if false these constraints are written in the
  /// form l_i <= h_i (x,y) <= u_i

  void writeMod (const std::string &fname, bool aux);

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

  /// fill vector with coefficients of objective function
  void fillObjCoeff (double *&);

  /// replace all occurrences of original variable with new aux given
  /// as argument
  void auxiliarize (exprAux *);
};


/// read best known solution from file given in argument
bool readOptimum (const std::string &, CouNumber *&, CouNumber &, CouenneProblem *);

#endif
