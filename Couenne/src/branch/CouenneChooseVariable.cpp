/*
 * Name:    CouenneChooseVariable.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneChooseVariable.hpp"
#include "CouenneProblem.hpp"


/// Default Constructor 
CouenneChooseVariable::CouenneChooseVariable (): 
  OsiChooseVariable (),
  problem_ (NULL) {}


/// Constructor from solver (so we can set up arrays etc)
CouenneChooseVariable::CouenneChooseVariable (const OsiSolverInterface *si,
					      CouenneProblem *p,
					      JnlstPtr jnlst):
  OsiChooseVariable (si),
  problem_ (p),
  jnlst_   (jnlst) {}


/// Copy constructor 
CouenneChooseVariable::CouenneChooseVariable (const CouenneChooseVariable &source):
  OsiChooseVariable (source),
  problem_ (source.problem_),
  jnlst_   (source.jnlst_) {}


/// Assignment operator 
CouenneChooseVariable & CouenneChooseVariable::operator= (const CouenneChooseVariable& rhs) {
  problem_ = rhs.problem_; 
  jnlst_   = rhs.jnlst_;
  return *this;
}


/// Sets up strong list and clears all if initialize is true.
/// Returns number of infeasibilities. 
/// If returns -1 then has worked out node is infeasible!
int CouenneChooseVariable::setupList (OsiBranchingInformation *info, bool initialize) {

  problem_ -> domain () -> push 
    (problem_ -> nVars (),
     info -> solution_, 
     info -> lower_, 
     info -> upper_);

  if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
    printf ("----------------- setup list\n");
    for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
      printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
	      info -> solution_ [i],
	      info -> lower_ [i],
	      info -> upper_ [i]);
  }

  // Make it stable, in OsiChooseVariable::setupList() numberObjects must be 0.
  int retval = (solver_ -> numberObjects ()) ? 
    OsiChooseVariable::setupList (info, initialize) : 0;

  problem_ -> domain () -> pop ();

  jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "----------------- setup list done\n");

  return retval;
}


/// choose object to branch based on earlier setup
// int CouenneChooseVariable::chooseVariable (OsiSolverInterface * solver, 
// 					   OsiBranchingInformation *info, 
// 					   bool fixVariables) {

//   // !!!should go -- just choose the first element
//   problem_ -> domain () -> push 
//     (problem_ -> nVars (),
//      info -> solution_, 
//      info -> lower_, 
//      info -> upper_);

//   int retval = OsiChooseVariable::chooseVariable (solver, info, fixVariables);

//   problem_ -> domain () -> pop ();

//   return retval;
// }


// Returns true if solution looks feasible against given objects
bool CouenneChooseVariable::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects) {
  double obj = solution [problem_ -> Obj (0) -> Body () -> Index ()];
  return problem_ -> checkNLP (solution, obj);
}


/// Add list of options to be read from file
void CouenneChooseVariable::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2 
    ("enable_sos",
     "Use Special Ordered Sets (SOS) as indicated in the MINLP model",
     "no",
     "no","",
     "yes","");

  roptions -> AddStringOption2
    ("branch_fbbt",
     "Apply bound tightening before branching",
     "yes",
     "no","",
     "yes","");

  roptions -> AddStringOption2
    ("branch_conv_cuts",
     "Apply convexification cuts before branching (for now only within strong branching)",
     "yes",
     "no","",
     "yes","");

  roptions -> AddStringOption6
    ("branch_pt_select",
     "Chooses branching point selection strategy",
     "mid-point",
     "lp-clamped", "LP point clamped in [k,1-k] of the bound intervals (k defined by lp_clamp)",
     "lp-central", "LP point if within [k,1-k] of the bound intervals, middle point otherwise" 
     "(k defined by branch_lp_clamp)",
     "balanced", "minimizes max distance from curve to convexification",
     "min-area", "minimizes total area of the two convexifications",
     "mid-point", "convex combination of current point and mid point",
     "no-branch", "do not branch, return null infeasibility; for testing purposes only",
     "");

  std::string br_ops [] = {"prod", "div", "exp", "log", "trig", 
			   "pow",  "negpow", "sqr", "cube", ""};

  for (int i=0; br_ops [i] != ""; i++) {

    char optname [40], optname2 [40], description [90];
    sprintf (optname,  "branch_pt_select_%s", br_ops [i].c_str ());
    sprintf (optname2, "branch_lp_clamp_%s",  br_ops [i].c_str ());
    sprintf (description, "Chooses branching point selection strategy for operator %s", 
	     br_ops [i].c_str ());

    roptions -> AddStringOption7
      (optname,
       description,
       "common",
       "common",    "use strategy defined for generic operators",
       "lp-clamped", "LP point clamped in [k,1-k] of the bound intervals "
       "(k defined by lp_clamp_${this operator})",
       "lp-central", "LP point if within [k,1-k] of the bound intervals, middle point otherwise" 
       "(k defined by branch_lp_clamp_${this operator})",
       "balanced",  "minimizes max distance from curve to convexification",
       "min-area",  "minimizes total area of the two convexifications",
       "mid-point", "convex combination of current point and mid point",
       "no-branch", "do not branch, return null infeasibility; for testing purposes only",
       "");

    roptions -> AddBoundedNumberOption
      (optname2,
       "Defines safe interval percentage [0,0.5] for using LP point as a branching point",
       0.,false,
       0.5,false,
       0.2,
       "Default value is 0.2.");
  }

  roptions -> AddBoundedNumberOption
    ("branch_midpoint_alpha",
     "Defines convex combination of mid point and current LP point: "
     "b = alpha x_lp + (1-alpha) (lb+ub)/2.",
     0.,false,
     1.,false,
     0.25,
     "Default value is 0.25.");

  roptions -> AddBoundedNumberOption
    ("branch_lp_clamp",
     "Defines safe interval percentage for using LP point as a branching point",
     0.,false,
     1.,false,
     0.2,
     "Default value is 0.2.");

  roptions -> AddLowerBoundedIntegerOption
    ("cont_var_priority",
     "Priority of continuous variable branching",
     1, 2000, "Default value is 2000 (it is 1000 for integers and 10 for SOS).");

  roptions -> AddStringOption2
    ("red_cost_branching",
     "Apply Reduced Cost Branching (instead of the Violation Transfer) -- MUST have vt_obj enabled",
     "no",
     "no", "Use Violation Transfer with $\\sum |\\pi_i a_{ij}|$",
     "yes","Use Reduced cost branching with $|\\sum \\pi_i a_{ij}|$");
}
