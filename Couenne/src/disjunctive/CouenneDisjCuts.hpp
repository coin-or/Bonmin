/*
 * Name:    CouenneDisjCuts.hpp
 * Author:  Pietro Belotti
 * Purpose: a generator of disjunctive cuts for MINLP problems
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_DISJUNCTIVE_CUTS_HPP
#define COUENNE_DISJUNCTIVE_CUTS_HPP

#include "BonRegisteredOptions.hpp"

#include "OsiSolverInterface.hpp"
#include "CglCutGenerator.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonBabSetupBase.hpp"
#include "BonBabInfos.hpp"
#include "OsiChooseVariable.hpp"
#include "CouenneTypes.hpp"
#include "CouenneJournalist.hpp"

class CouenneProblem;
class CouenneSolverInterface;

enum {COUENNE_INFEASIBLE, COUENNE_TIGHTENED, COUENNE_FEASIBLE};

/// Cut Generator for linear convexifications

class CouenneDisjCuts: public CglCutGenerator {

 protected:

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
  /// Couenne pass of each b&b node)
  Bonmin::OsiTMINLPInterface *minlp_;

  /// Branching scheme (if strong, we can use SB candidates)
  OsiChooseVariable *branchingMethod_;

  /// Is branchMethod_ referred to a strong branching scheme?
  bool isBranchingStrong_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// Number of disjunction to consider at each separation 
  mutable int numDisjunctions_;

  /// Initial percentage of objects to use for generating cuts, in [0,1]
  double initDisjPercentage_;

  /// Depth of the BB tree where start decreasing number of objects
  int depthLevelling_;

  /// Depth of the BB tree where stop separation
  int depthStopSeparate_;

 public:

  /// constructor
  CouenneDisjCuts (Bonmin::OsiTMINLPInterface *minlp = NULL,
		   Bonmin::BabSetupBase *base = NULL,
		   CouenneProblem *problem = NULL,
		   OsiChooseVariable *bcv = NULL,
		   bool is_strong = false,
		   JnlstPtr journalist = NULL,
		   const Ipopt::SmartPtr<Ipopt::OptionsList> options = NULL);

  /// copy constructor
  CouenneDisjCuts (const CouenneDisjCuts &);

  /// destructor
  ~CouenneDisjCuts ();

  /// clone method (necessary for the abstract CglCutGenerator class)
  CouenneDisjCuts *clone () const
  {return new CouenneDisjCuts (*this);}

  /// return pointer to symbolic problem
  inline CouenneProblem *Problem () const
    {return problem_;}

  /// the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  /// Provide Journalist
  inline ConstJnlstPtr Jnlst() const 
  {return ConstPtr (jnlst_);}

  /// get all disjunctions
  int getDisjunctions (std::vector <std::pair <OsiCuts *, OsiCuts *> > &disjunctions, 
		       OsiSolverInterface &si, 
		       OsiCuts &cs, 
		       const CglTreeInfo &info) const;

  /// separate couenne cuts on both sides of single disjunction
  int separateWithDisjunction (std::pair <OsiCuts *, OsiCuts *> &disj, 
			      OsiSolverInterface &si, 
			      OsiCuts &cs, 
			      const CglTreeInfo &info) const;

  /// generate one disjunctive cut from one CGLP
  int generateDisjCut (std::pair <OsiCuts *, OsiCuts *> &disj, 
		       OsiSolverInterface &si, 
		       OsiCuts &cs, 
		       const CglTreeInfo &info) const;

  /// check if (column!) cuts compatible with solver interface
  int checkDisjSide (OsiSolverInterface &si, OsiCuts *cuts) const;

  /// compute smallest box containing both left and right boxes.
  int getBoxUnion (OsiSolverInterface &si, 
		   OsiCuts *left, OsiCuts *right, 
		   CoinPackedVector &lower, CoinPackedVector &upper) const;

protected:

  /// create single osicolcut disjunction
  OsiCuts *getSingleDisjunction (OsiSolverInterface &si) const;
};

#endif
