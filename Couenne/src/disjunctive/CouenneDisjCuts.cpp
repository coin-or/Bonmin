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
				  CouenneCutGenerator *cg,
				  OsiChooseVariable *bcv,
				  bool is_strong,
				  JnlstPtr journalist,
				  const Ipopt::SmartPtr<Ipopt::OptionsList> options):
  couenneCG_          (cg),
  nrootcuts_          (-1), // to indicate first iteration not done yet
  ntotalcuts_         (0),
  septime_            (0.),
  objValue_           (-COIN_DBL_MAX),
  minlp_              (minlp),
  branchingMethod_    (bcv),
  isBranchingStrong_  (is_strong),
  jnlst_              (journalist),
  activeRows_         (false),
  activeCols_         (false),
  addPreviousCut_     (false),
  cpuTime_            (-1.) {

  options -> GetNumericValue ("time_limit", cpuTime_,  "couenne.");

  options -> GetNumericValue ("disj_init_perc",   initDisjPercentage_,  "couenne.");
  options -> GetIntegerValue ("disj_init_number", initDisjNumber_,      "couenne.");
  options -> GetIntegerValue ("disj_depth_level", depthLevelling_,      "couenne.");
  options -> GetIntegerValue ("disj_depth_stop",  depthStopSeparate_,   "couenne.");

  std::string s;
  options -> GetStringValue ("disj_active_rows", s, "couenne."); activeRows_     = (s == "yes");
  options -> GetStringValue ("disj_active_cols", s, "couenne."); activeCols_     = (s == "yes");
  options -> GetStringValue ("disj_cumulative",  s, "couenne."); addPreviousCut_ = (s == "yes");
}


/// copy constructor
CouenneDisjCuts::CouenneDisjCuts (const CouenneDisjCuts &src):
  couenneCG_          (src.couenneCG_),
  nrootcuts_          (src.nrootcuts_),
  ntotalcuts_         (src.ntotalcuts_),
  septime_            (src.septime_),
  objValue_           (src.objValue_),
  minlp_              (src.minlp_),
  branchingMethod_    (src.branchingMethod_),
  isBranchingStrong_  (src.isBranchingStrong_),
  jnlst_              (src.jnlst_),
  initDisjPercentage_ (src.initDisjPercentage_),
  initDisjNumber_     (src.initDisjNumber_),
  depthLevelling_     (src.depthLevelling_),
  depthStopSeparate_  (src.depthStopSeparate_),
  activeRows_         (src.activeRows_),
  activeCols_         (src.activeCols_),
  addPreviousCut_     (src.addPreviousCut_),
  cpuTime_            (src.cpuTime_) {}


/// Add list of options to be read from file
void CouenneDisjCuts::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddLowerBoundedIntegerOption
    ("minlp_disj_cuts",
     "The frequency (in terms of nodes) at which Couenne disjunctive cuts are generated.",
     0, 0,
     "A frequency of 0 (default) means these cuts are never generated.");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_init_number",
     "Maximum number of disjunction to consider at each iteration (-1 for unlimited, default 10).",
     -1, 10, "");

  roptions -> AddBoundedNumberOption
    ("disj_init_perc",
     "The maximum fraction of OsiObjects to consider for generating disjunctions",
     0., false,
     1., false,
     0.5, "");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_level",
     "Depth of the BB tree when to start decreasing the number of objects that generate disjunctions.",
     -1, 5, "");

  roptions -> AddLowerBoundedIntegerOption
    ("disj_depth_stop",
     "Depth of the BB tree where separation is stopped.",
     -1, 20, "");

  roptions -> AddStringOption2
    ("disj_active_rows",
     "Only include active rows in the CGLP.",
     "no", 
     "yes", "",
     "no", "");

  roptions -> AddStringOption2
    ("disj_active_cols",
     "Only include active cols in the CGLP.",
     "no", 
     "yes", "",
     "no", "");

  roptions -> AddStringOption2
    ("disj_cumulative",
     "Add previous disj. cut to current CGLP.",
     "no", 
     "yes", "",
     "no", "");
}
