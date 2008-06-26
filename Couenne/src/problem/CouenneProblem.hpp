/*
 * Name:    CouenneProblem.hpp
 * Author:  Pietro Belotti
 * Purpose: define the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PROBLEM_HPP
#define COUENNE_PROBLEM_HPP

#include <vector>

#include "OsiRowCut.hpp"
#include "BonBabInfos.hpp"

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprAux.hpp"
#include "CouenneProblemElem.hpp"
#include "CouenneObject.hpp"
#include "CouenneJournalist.hpp"
#include "domain.hpp"

struct ASL;
struct expr;

class DepGraph;
class CouenneCutGenerator;
class CouenneSolverInterface;
class quadElem;
class LinMap;
class QuadMap;

// default tolerance for checking feasibility (and integrality) of NLP solutions
const CouNumber feas_tolerance_default = 1e-5;


/** Class for MINLP problems with symbolic information
 *
 *  It is read from an AMPL .nl file and contains variables, AMPL's
 *  "defined variables" (aka common expressions), objective(s), and
 *  constraints in the form of expression's. Changes throughout the
 *  program occur in standardization.
 */

class CouenneProblem {

  /// Class for storing a global cutoff for a CouenneProblem and all
  /// its clones
  class GlobalCutOff {
  private:
    GlobalCutOff(const GlobalCutOff&);
    double cutoff_;
  public:
    GlobalCutOff ()         : cutoff_ (COIN_DBL_MAX) {}
    GlobalCutOff (double c) : cutoff_ (c) {}
    ~GlobalCutOff() {}
    inline void setCutOff (double cutoff) {cutoff_ = cutoff;}
    inline double getCutOff() const {return cutoff_;}
  };

  /// structure to record fixed, non-fixed, and continuous variables
  enum fixType {UNFIXED, FIXED, CONTINUOUS};

 protected:

  /// problem name
  std::string problemName_;

  std::vector <exprVar           *> variables_;   ///< Variables (original, auxiliary, and defined)
  std::vector <CouenneObjective  *> objectives_;  ///< Objectives
  std::vector <CouenneConstraint *> constraints_; ///< Constraints

  /// AMPL's common expressions (read from AMPL through structures cexps and cexps1)
  std::vector <expression *> commonexprs_; 

  mutable Domain domain_; ///< current point and bounds;

  /// Expression map for comparison in standardization and to count
  /// occurrences of an auxiliary
  std::set <exprAux *, compExpr> *auxSet_;

  /// Number of elements in the x_, lb_, ub_ arrays
  mutable int curnvars_;

  /// Number of discrete variables
  int nIntVars_;

  /// Best solution known to be loaded from file -- for testing purposes
  CouNumber *optimum_;

  /// Best known objective function
  CouNumber bestObj_;

  /// Indices of variables appearing in products (used for SDP cuts)
  int *quadIndex_;

  /// Variables that have commuted to auxiliary
  bool *commuted_;

  /// numbering of variables. No variable xi with associated pi(i)
  /// greater than pi(j) should be evaluated before variable xj
  int *numbering_;

  /// Number of "defined variables" (aka "common expressions")
  int ndefined_;

  /// Dependence (acyclic) graph: shows dependence of all auxiliary
  /// variables on one another and on original variables. Used to
  /// create a numbering of all variables for evaluation and bound
  /// tightening (reverse order for implied bounds)
  DepGraph *graph_;

  /// Number of original variables
  int nOrig_;

  /// Number of original constraints (disregarding those that turned
  /// into auxiliary variable definition)
  int nOrigCons_;

  /// Pointer to a global cutoff object
  mutable GlobalCutOff* pcutoff_;

  /// flag indicating if this class is creator of global cutoff object
  mutable bool created_pcutoff_;

  bool doFBBT_;  ///< do Feasibility-based bound tightening
  bool doOBBT_;  ///< do Optimality-based  bound tightening
  bool doABT_;   ///< do Aggressive        bound tightening

  int logObbtLev_;   ///< frequency of Optimality-based bound tightening
  int logAbtLev_;    ///< frequency of Aggressive       bound tightening

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// window around known optimum (for testing purposes)
  CouNumber opt_window_;

  /// Use quadratic expressions?
  bool useQuadratic_;

  /// feasibility tolerance (to be used in checkNLP)
  CouNumber feas_tolerance_;

  /// inverse dependence structure: for each variable x give set of
  /// auxiliary variables (or better, their indices) whose expression
  /// depends on x
  std::vector <std::set <int> > dependence_;

  /// vector of pointer to CouenneObjects. Used by CouenneVarObjects
  /// when finding all objects related to (having as argument) a
  /// single variable
  std::vector <CouenneObject> objects_;

  /// each element is true if variable is integer and, if auxiliary,
  /// depends on no integer
  mutable int *integerRank_;

  /// numberInRank_ [i] is the number of integer variables in rank i
  mutable std::vector <int> numberInRank_;

 public:

  CouenneProblem  (const ASL * = NULL,
		   Bonmin::BabSetupBase *base = NULL,
		   JnlstPtr jnlst = NULL);  ///< Constructor
  CouenneProblem  (const CouenneProblem &); ///< Copy constructor
  ~CouenneProblem ();                       ///< Destructor

  /// Clone method (for use within CouenneCutGenerator::clone)
  CouenneProblem *clone () const
  {return new CouenneProblem (*this);}

  int nObjs     () const {return objectives_.   size ();} ///< Get number of objectives
  int nCons     () const {return constraints_.  size ();} ///< Get number of constraints
  int nOrigCons () const {return nOrigCons_;}             ///< Get number of original constraints

  inline int nOrig    () const {return nOrig_;}              ///< Number of original variables
  inline int nIntVars () const {return nIntVars_;}           ///< Number of original integer variables
  inline int nVars    () const {return variables_. size ();} ///< Total number of variables

  /// get evaluation order index 
  inline int evalOrder (int i) const
  {return numbering_ [i];}

  /// get evaluation order vector (numbering_)
  inline int *evalVector ()
  {return numbering_;}

  // get elements from vectors
  inline CouenneConstraint *Con (int i) const {return constraints_ [i];} ///< i-th constraint
  inline CouenneObjective  *Obj (int i) const {return objectives_  [i];} ///< i-th objective

  /// Return pointer to i-th variable
  inline exprVar *Var   (int i) const 
  {return variables_ [i];}

  /// Return vector of variables (symbolic representation)
  inline std::vector <exprVar *> &Variables () 
  {return variables_;}

  /// Return pointer to set for comparisons
  inline std::set <exprAux *, compExpr> *AuxSet () 
  {return auxSet_;}

  /// Return pointer to dependence graph
  inline DepGraph *getDepGraph () 
  {return graph_;}

  /// return current point & bounds
  inline Domain *domain () const
  {return &domain_;}

  // Get and set current variable and bounds
  inline CouNumber   &X     (int i) const {return domain_.x   (i);} ///< \f$x_i\f$
  inline CouNumber   &Lb    (int i) const {return domain_.lb  (i);} ///< lower bound on \f$x_i\f$
  inline CouNumber   &Ub    (int i) const {return domain_.ub  (i);} ///< upper bound on\f$x_i\f$

  // get and set current variable and bounds
  inline CouNumber  *X     () const {return domain_.x  ();} ///< Return vector of variables
  inline CouNumber  *Lb    () const {return domain_.lb ();} ///< Return vector of lower bounds
  inline CouNumber  *Ub    () const {return domain_.ub ();} ///< Return vector of upper bounds

  // get optimal solution and objective value
  CouNumber  *bestSol () const {return optimum_;} ///< Best known solution (read from file)
  CouNumber   bestObj () const {return bestObj_;} ///< Objective of best known solution

  /// Get vector of commuted variables
  bool *&Commuted () 
  {return commuted_;}

  /// Add (non linear) objective function
  void addObjective     (expression *, const std::string &);

  // Add (non linear) "=", ">=", "<=", and range constraints
  void addEQConstraint  (expression *, expression *); ///< Add equality constraint \f$ h(x) = b\f$
  void addGEConstraint  (expression *, expression *); ///< Add \f$\ge\f$ constraint, \f$h(x)\ge b\f$
  void addLEConstraint  (expression *, expression *); ///< Add \f$\le\f$ constraint, \f$h(x)\le b\f$
  void addRNGConstraint (expression *, expression *, 
			 expression *);               ///< Add range constraint, \f$a\le h(x)\le b\f$

  /// Add original variable.
  ///
  /// @param isint if true, this variable is integer, otherwise it is
  /// continuous
  expression *addVariable (bool isint = false, Domain *d = NULL);

  /// Add auxiliary variable and associate it with expression given as
  /// argument (used in standardization)
  exprAux *addAuxiliary (expression *);

  /// Break problem's nonlinear constraints in simple expressions to be
  /// convexified later
  void standardize ();

  /// Display current representation of problem: objective, linear and
  /// nonlinear constraints, and auxiliary variables.
  void print (std::ostream & = std::cout);

  /// Read problem from .nl file using the Ampl Solver Library (ASL)
  int readnl (const struct ASL *);

  /// Generate a Couenne expression from an ASL expression
  expression *nl2e (struct expr *, const ASL *asl);

  // bound tightening parameters
  bool doFBBT () const {return doFBBT_;} ///< shall we do Feasibility Based Bound Tightening?
  bool doOBBT () const {return doOBBT_;} ///< shall we do Optimality  Based Bound Tightening?
  bool doABT  () const {return doABT_;}  ///< shall we do Aggressive        Bound Tightening?

  int  logObbtLev () const {return logObbtLev_;} ///< How often shall we do OBBT?
  int  logAbtLev  () const {return logAbtLev_;}  ///< How often shall we do ABT?

  /// Write nonlinear problem to a .mod file (with lots of defined
  /// variables)
  /// 
  /// @param fname Name of the .mod file to be written
  ///
  /// @param aux controls the use of auxiliaries. If true, a problem
  /// is written with auxiliary variables written with their
  /// associated expression, i.e. \f$w_i = h_i(x,y,w)\f$ and bounds
  /// \f$l_i \le w_i \le u_i\f$, while if false these constraints are
  /// written in the form \f$l_i \le h_i (x,y) \le u_i\f$.
  ///
  /// Note: if used before standardization, writes original AMPL formulation
  void writeAMPL (const std::string &fname, bool aux);

  /// Write nonlinear problem to a .gms file
  /// 
  /// @param fname Name of the .gams file to be written.
  void writeGAMS (const std::string &fname);

  /// Initialize auxiliary variables and their bounds from original
  /// variables
  //void initAuxs (const CouNumber *, const CouNumber *, const CouNumber *);
  void initAuxs () const;

  /// Get auxiliary variables from original variables
  void getAuxs (CouNumber *) const;

  /// tighten bounds using propagation, implied bounds and reduced costs
  bool boundTightening (t_chg_bounds *, 
			Bonmin::BabInfo * = NULL) const;

  /// core of the bound tightening procedure
  bool btCore (t_chg_bounds *chg_bds) const;

  /// Optimality Based Bound Tightening
  int obbt (const CouenneCutGenerator *cg,
	    const OsiSolverInterface &csi,
	    OsiCuts &cs,
	    const CglTreeInfo &info,
	    Bonmin::BabInfo * babInfo,
	    t_chg_bounds *chg_bds);

  /// aggressive bound tightening. Fake bounds in order to cut
  /// portions of the solution space by fathoming on
  /// bounds/infeasibility
  bool aggressiveBT (Bonmin::OsiTMINLPInterface *nlp,
		     t_chg_bounds *, 
		     Bonmin::BabInfo * = NULL) const;

  /// procedure to strengthen variable bounds. Return false if problem
  /// turns out to be infeasible with given bounds, true otherwise.
  int redCostBT (const OsiSolverInterface *psi,
		 t_chg_bounds *chg_bds) const;

  /// "Forward" bound tightening, that is, propagate bound of variable
  /// \f$x\f$ in an expression \f$w = f(x)\f$ to the bounds of \f$w\f$.
  int tightenBounds (t_chg_bounds *) const;

  /// "Backward" bound tightening, aka implied bounds. 
  int impliedBounds (t_chg_bounds *) const;

  /// Look for quadratic terms to be used with SDP cuts
  void fillQuadIndices ();

  /// Fill vector with coefficients of objective function
  void fillObjCoeff (double *&);

  /// Replace all occurrences of original variable with new aux given
  /// as argument
  void auxiliarize (exprVar *, exprVar * = NULL);

  /// Set cutoff
  void setCutOff (CouNumber cutoff) const;

  /// Set cutoff
  CouNumber getCutOff () const
  {return pcutoff_ -> getCutOff ();}

  /// Make cutoff known to the problem
  void installCutOff ();

  /// Provide Journalist
  ConstJnlstPtr Jnlst() const 
  {return ConstPtr(jnlst_);}

  /// Check if solution is MINLP feasible
  bool checkNLP (const double *solution, const double obj) const;

  /// generate integer NLP point Y starting from fractional solution
  /// using bound tightening
  int getIntegerCandidate (const double *xFrac, double *xInt, double *lb, double *ub) const;

  /// Read best known solution from file given in argument
  bool readOptimum (std::string *fname = NULL);

  /// Read cutoff value (for testing purposes)
  void readCutoff (const std::string &fname);

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  /// standardization of linear exprOp's
  exprAux *linStandardize (bool addAux, 
			   CouNumber c0, 
			   LinMap  &lmap,
			   QuadMap &qmap);

  /// split a constraint w - f(x) = c into w's index (it is returned)
  /// and rest = f(x) + c
  int splitAux (CouNumber, expression *, expression *&, bool *);

  /// translates pair (indices, coefficients) into vector with pointers to variables
  void indcoe2vector (int *indexL,
		      CouNumber *coeff,
		      std::vector <std::pair <exprVar *, CouNumber> > &lcoeff);

  /// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
  void indcoe2vector (int *indexI,
		      int *indexJ,
		      CouNumber *coeff,
		      std::vector <quadElem> &qcoeff);

  /// given (expression *) element of sum, returns (coe,ind0,ind1)
  /// depending on element:
  ///
  /// 1) a * x_i ^ 2   ---> (a,i,?)   return COU_EXPRPOW
  /// 2) a * x_i       ---> (a,i,?)   return COU_EXPRVAR
  /// 3) a * x_i * x_j ---> (a,i,j)   return COU_EXPRMUL
  /// 4) a             ---> (a,?,?)   return COU_EXPRCONST
  ///
  /// x_i and/or x_j may come from standardizing other (linear or
  /// quadratic operator) sub-expressions
  void decomposeTerm (expression *term,
		      CouNumber initCoe,
		      CouNumber &c0,
		      LinMap  &lmap,
		      QuadMap &qmap);

  /// return problem name
  const std::string &problemName () const
  {return problemName_;}

  /// return inverse dependence structure
  const std::vector <std::set <int> > &Dependence () const
  {return dependence_;}

  /// return object vector
  const std::vector <CouenneObject> &Objects () const
  {return objects_;}

  /// find SOS constraints in problem
  int findSOS (OsiSolverInterface *solver, OsiObject ** objects);

protected:

  /// single fake tightening. Return
  ///
  /// -1   if infeasible
  ///  0   if no improvement
  /// +1   if improved
  int fake_tighten (char direction,  ///< 0: left, 1: right
		    int index,       ///< index of the variable tested
		    const double *X, ///< point round which tightening is done
		    CouNumber *olb,  ///< cur. lower bound
		    CouNumber *oub,  ///< cur. upper bound
		    t_chg_bounds *chg_bds,
		    t_chg_bounds *f_chg) const;

  /// Optimality Based Bound Tightening -- inner loop
  int obbtInner (CouenneSolverInterface *, 
		 OsiCuts &,
		 t_chg_bounds *, 
		 Bonmin::BabInfo *) const;

  int obbt_iter (CouenneSolverInterface *csi, 
		 t_chg_bounds *chg_bds, 
		 const CoinWarmStart *warmstart, 
		 Bonmin::BabInfo *babInfo,
		 double *objcoe,
		 int sense, 
		 int index) const;

  int call_iter (CouenneSolverInterface *csi, 
		 t_chg_bounds *chg_bds, 
		 const CoinWarmStart *warmstart, 
		 Bonmin::BabInfo *babInfo,
		 double *objcoe,
		 enum nodeType type,
		 int sense) const;

  /// analyze sparsity of potential exprQuad/exprGroup and change
  /// linear/quadratic maps accordingly, if necessary by adding new
  /// auxiliary variables and including them in the linear map
  void analyzeSparsity (CouNumber, 
			LinMap &,
			QuadMap &);

  /// re-organizes multiplication and stores indices (and exponents) of
  /// its variables
  void flattenMul (expression *mul, 
		   CouNumber &coe, 
		   std::map <int, CouNumber> &indices);

  /// clear all spurious variables pointers not referring to the variables_ vector
  void realign ();

  /// fill dependence_ structure
  void fillDependence (Bonmin::BabSetupBase *base);

  /// fill freeIntegers_ array
  void fillIntegerRank () const;

  /// Test fixing of an integer variable (used in getIntegerCandidate())
  int testIntFix (int index, 
		  CouNumber xFrac, 
		  enum fixType *fixed,
		  CouNumber *xInt,
		  CouNumber *dualL, CouNumber *dualR,
		  CouNumber *olb,   CouNumber *oub,
		  bool patient) const;
};

#endif
