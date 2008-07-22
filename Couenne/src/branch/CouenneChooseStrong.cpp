/*
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneBranchingObject.hpp"

namespace Bonmin {

  /// constructor
  CouenneChooseStrong::CouenneChooseStrong (BabSetupBase &b, CouenneProblem* p, JnlstPtr jnlst) :

    BonChooseVariable (b, b.continuousSolver()),
    problem_          (p),
    jnlst_            (jnlst) {

    std::string s;
    b.options () -> GetStringValue ("pseudocost_mult_lp", s, "couenne.");
    pseudoUpdateLP_ = (s == "yes");      
  }

  /// copy constructor
  CouenneChooseStrong::CouenneChooseStrong (const CouenneChooseStrong& rhs) :
    BonChooseVariable (rhs),
    problem_          (rhs.problem_),
    pseudoUpdateLP_   (rhs.pseudoUpdateLP_),
    jnlst_            (rhs.jnlst_)
  {}

  /// destructor
  CouenneChooseStrong::~CouenneChooseStrong()
  {}

  /// cloning method
  OsiChooseVariable *
  CouenneChooseStrong::clone() const
  {
    return new CouenneChooseStrong(*this);
  }

  /// assignment operator
  CouenneChooseStrong&
  CouenneChooseStrong::operator=(const CouenneChooseStrong & rhs)
  {
    if (this != &rhs) {
      BonChooseVariable::operator=(rhs);
      problem_ = rhs.problem_;
    }
    return *this;
  }


  /// choose object to branch based on earlier setup
  int CouenneChooseStrong::chooseVariable (OsiSolverInterface * solver,
					   OsiBranchingInformation *info,
					   bool fixVariables) {
    problem_ -> domain () -> push
      (problem_ -> nVars (),
       info -> solution_,
       info -> lower_,
       info -> upper_);

    int retval = BonChooseVariable::chooseVariable (solver, info, fixVariables);

    problem_ -> domain () -> pop ();

    return retval;
  }


  // Sets up strong list and clears all if initialize is true.
  // Returns number of infeasibilities.
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {

    initialize = true; // to avoid failed assert in BonChooseVariable::setupList()

    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       info -> solution_, 
       info -> lower_, 
       info -> upper_); // have to alloc+copy

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("----------------- (strong) setup list\n");
      for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
	printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
		info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
    }

    // call Bonmin's setuplist
    int retval = BonChooseVariable::setupList (info, initialize);

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, 
		      "----------------- (strong) setup list done - %d infeasibilities\n", retval);

    problem_ -> domain () -> pop ();
    return retval;
  }


  /// Add list of options to be read from file ////////////////////////////////////////
  void CouenneChooseStrong::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

    roptions -> AddStringOption6
      ("pseudocost_mult",
       "Multipliers of pseudocosts for estimating and update estimation of bound",
       "interval_br_rev",

       "infeasibility", "infeasibility returned by object",

       "projectDist",   "distance between current LP point and resulting branches' LP points",

       "interval_lp",   "width of the interval between bound and current lp point",
       "interval_lp_rev",   "similar to interval_lp, reversed",

       "interval_br",   "width of the interval between bound and branching point",
       "interval_br_rev",   "similar to interval_br, reversed");

    roptions -> AddStringOption2
      ("pseudocost_mult_lp",
       "Use distance between LP points to update multipliers of pseudocosts "  
       "after simulating branching",
       "no",
       "yes", "",
       "no",  "");
  }


  // Returns true if solution looks feasible against given objects
  bool CouenneChooseStrong::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects)

  {return problem_ -> checkNLP (solution, solution [problem_ -> Obj (0) -> Body () -> Index ()]);}
}
