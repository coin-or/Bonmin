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

 public:

  /// constructor
  CouenneDisjCuts (Bonmin::OsiTMINLPInterface *minlp = NULL,
		   Bonmin::BabSetupBase *base = NULL,
		   CouenneProblem *problem = NULL,
		   OsiChooseVariable *bcv = NULL,
		   bool is_strong = false,
		   JnlstPtr journalist = NULL);

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
};

#endif
