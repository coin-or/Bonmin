/*
 * Name:    CouenneDisjCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for the disjunctive cuts
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"


/// constructor
CouenneDisjCuts::CouenneDisjCuts (Bonmin::OsiTMINLPInterface *minlp,
				  Bonmin::BabSetupBase *base,
				  CouenneProblem *problem,
				  OsiChooseVariable *bcv,
				  bool is_strong,
				  JnlstPtr journalist,
				  const Ipopt::SmartPtr<Ipopt::OptionsList> options):

  problem_            (problem),
  nrootcuts_          (-1), // to indicate first iteration not done yet
  ntotalcuts_         (0),
  septime_            (0.),
  objValue_           (-COIN_DBL_MAX),
  minlp_              (minlp),
  branchingMethod_    (bcv),
  isBranchingStrong_  (is_strong),
  jnlst_              (journalist) {

  options -> GetNumericValue ("disj_init_perc",   initDisjPercentage_,  "couenne.");
  options -> GetIntegerValue ("disj_depth_level", depthLevelling_,      "couenne.");
  options -> GetIntegerValue ("disj_depth_stop",  depthStopSeparate_,   "couenne.");
}


/// copy constructor
CouenneDisjCuts::CouenneDisjCuts (const CouenneDisjCuts &src):
  problem_            (src.problem_),
  nrootcuts_          (src.nrootcuts_),
  ntotalcuts_         (src.ntotalcuts_),
  septime_            (src.septime_),
  objValue_           (src.objValue_),
  minlp_              (src.minlp_),
  branchingMethod_    (src.branchingMethod_),
  isBranchingStrong_  (src.isBranchingStrong_),
  jnlst_              (src.jnlst_),
  initDisjPercentage_ (src.initDisjPercentage_),
  depthLevelling_     (src.depthLevelling_),
  depthStopSeparate_  (src.depthStopSeparate_) {}


/// destructor
CouenneDisjCuts::~CouenneDisjCuts () {}


/// Add list of options to be read from file
void CouenneDisjCuts::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("minlp_disj_cuts",
     "The frequency (in terms of nodes) at which Couenne disjunctive cuts are generated.",
     0, 0,
     "A frequency of 0 (default) means these cuts are never generated.");

  roptions -> AddBoundedNumberOption
    ("disj_init_perc",
     "The maximum fraction of OsiObjects to consider for generating disjunctions",
     0., false,
     1., false,
     0.5, "");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_level",
     "Depth of the BB tree when we start decreasing the number of objects to generate disjunctions.",
     0, 5, "");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_stop",
     "Depth of the BB tree when we stop separation.",
     -1, 20, "");
}
