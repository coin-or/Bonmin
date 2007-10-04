/*
 * Name: CouenneCglCutGenerator.hpp
 * Author: Pietro Belotti
 * Purpose: define a class of pools of convexification OsiRowCuts
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#ifndef COUENNE_CUT_GENERATOR_HPP
#define COUENNE_CUT_GENERATOR_HPP

#include <BonOaDecBase.hpp>
#include <OsiRowCut.hpp>
#include <CouenneTypes.hpp>
#include "BonAuxInfos.hpp"

#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>

#include <BonCbc.hpp>

class CouenneProblem;

struct ASL;

/// Cut Generator for linear convexifications

class CouenneCutGenerator: public Bonmin::OaDecompositionBase {

 protected:

  /// has generateCuts been called before?
  mutable bool firstcall_;

  /// should we add the violated cuts only (true), or all of them (false)?
  mutable bool addviolated_;

  /// what kind of sampling should be performed?
  enum conv_type convtype_;

  /// how many cuts should be added for each function?
  int nSamples_;

  /// pointer to symbolic repr. of constraint, variables, and bounds
  CouenneProblem *problem_;

  /// number of cuts generated at the first call
  mutable int nrootcuts_;

  /// total number of cuts generated 
  mutable int ntotalcuts_;

  /// separation time (includes generation of problem)
  mutable double septime_;

  /// Record obj value at final point of CouenneConv.
  mutable double objValue_;

  /// nonlinear solver interface as used within Bonmin (used at first
  /// Couenne pass of each b&b node
  Bonmin::OsiTMINLPInterface *nlp_;

  /// pointer to the Bab object (used to retrieve the current primal
  /// bound through bestObj())
  Bonmin::Bab *BabPtr_;

  /// signal infeasibility of current node (found through bound tightening)
  mutable bool infeasNode_;

 public:

  /// constructor
  CouenneCutGenerator  (Bonmin::OsiTMINLPInterface * = NULL,
			const struct ASL * = NULL, 
			bool = false, 
			enum conv_type = UNIFORM_GRID, 
			int = 2);

  /// copy constructor
  CouenneCutGenerator  (const CouenneCutGenerator &);

  /// destructor
  ~CouenneCutGenerator ();

  /// clone method (necessary for the abstract CglCutGenerator class)
  CouenneCutGenerator *clone () const
  {return new CouenneCutGenerator (*this);}

  /// return pointer to symbolic problem
  CouenneProblem *Problem () const
    {return problem_;}

  /// total number of variables (original + auxiliary)
  const int getnvars () const;

  /// has generateCuts been called yet?
  bool isFirst () const 
    {return firstcall_;}

  /// should we add the violated cuts only (true), or all of them (false)?
  bool addViolated () const
    {return addviolated_;}

  /// get convexification type (see CouenneTypes.h)
  enum conv_type ConvType () const
    {return convtype_;}

  /// get number of convexification samples
  int nSamples () const
    {return nSamples_;}

  /// the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  /// create cut and check violation. Insert and return status
  int createCut (OsiCuts &, // cutset to insert
		 CouNumber, // lb
		 CouNumber, // ub
		            // index, coeff  (index -1: "don't care") 
		 int,    CouNumber,    // of first  term
		 int=-1, CouNumber=0., // of second term 
		 int=-1, CouNumber=0., // of third  term
		 bool = false) const;  // is it a global cut? No, by default

  /// create cut and check violation. Other version with only one bound
  int createCut (OsiCuts &, // cutset to insert
		 CouNumber, // rhs
		 int,       // sign: -1: <=, 0: =, +1: >=
		            // index, coeff  (index -1: "don't care") 
		 int,    CouNumber,    // of first  term
		 int=-1, CouNumber=0., // of second term 
		 int=-1, CouNumber=0., // of third  term
		 bool = false) const;  // is it a global cut? No, by default

  /// Add general linear envelope to convex function, given its
  /// variables' indices, the (univariate) function and its first
  /// derivative
  void addEnvelope (OsiCuts &,
		    int,
		    unary_function, unary_function, 
		    int, int, 
		    CouNumber, CouNumber, CouNumber,
		    t_chg_bounds * = NULL,
		    bool = false) const;

  /// Add half-plane through (x1,y1) and (x2,y2) -- resp. 4th, 5th,
  /// 6th, and 7th argument
  int addSegment (OsiCuts &, int, int,
		  CouNumber, CouNumber, 
		  CouNumber, CouNumber, int) const;

  /// add tangent at given poing (x,w) with given slope
  int addTangent (OsiCuts &, int, int, 
		  CouNumber, CouNumber, 
		  CouNumber, int) const; 

  /// (for compatibility with base class)
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  virtual double performOa (OsiCuts & cs,           solverManip &nlpManip, 
			    solverManip &lpManip,   SubMipSolver *& subMip, 
			    OsiBabSolver * babInfo, double &cutoff) const
    {throw -1;}

  /// (for compatibility with base class)
  /// virtual method to decide if local search is performed
  virtual bool doLocalSearch () const {return 0;}

  /// Method to set the Bab pointer
  void setBabPtr (Bonmin::Bab *p)
    {BabPtr_ = p;}

  /// Get statistics
  void getStats (int &nrc, int &ntc, double &st) {
    nrc = nrootcuts_;
    ntc = ntotalcuts_;
    st  = septime_;
  }

  /// Allow to get and set the infeasNode_ flag (used only in generateCuts())
  bool &infeasNode () const
    {return infeasNode_;}

  /// generate OsiRowCuts for current convexification
  void genRowCuts (const OsiSolverInterface &, OsiCuts &cs, 
		   int, int *, const CglTreeInfo &, 
		   t_chg_bounds * = NULL, bool = false) const;

  /// generate OsiColCuts for improved (implied and propagated) bounds
  void genColCuts (const OsiSolverInterface &, OsiCuts &, int, int *) const;

  /// tighten bounds using propagation, implied bounds and reduced costs
  bool boundTightening (const OsiSolverInterface *, OsiCuts &, 
			t_chg_bounds *, Bonmin::BabInfo * = NULL) const;

  /// Optimality Based Bound Tightening
  int obbt (OsiSolverInterface *, OsiCuts &, t_chg_bounds *, Bonmin::BabInfo *) const;
};

#endif
