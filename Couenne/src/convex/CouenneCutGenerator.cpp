/*
 * Name: CouenneCglCutGenerator.cpp
 * Author: Pietro Belotti
 * Purpose: define a class of convexification procedures 
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */


#include "OsiRowCut.hpp"
#include "BonOaDecBase.hpp"
#include "CglCutGenerator.hpp"
#include "BonBabInfos.hpp"

#include "CouennePrecisions.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"


/// constructor
CouenneCutGenerator::CouenneCutGenerator (Bonmin::OsiTMINLPInterface *nlp,
					  Bonmin::BabSetupBase *base,
					  const struct ASL *asl, 
					  JnlstPtr jnlst):

  OaDecompositionBase (nlp, NULL, NULL, 0,0,0),

  firstcall_      (true),
  problem_        (NULL),
  nrootcuts_      (0),
  ntotalcuts_     (0),
  septime_        (0),
  objValue_       (- DBL_MAX),
  nlp_            (nlp),
  BabPtr_         (NULL),
  infeasNode_     (false),
  jnlst_          (jnlst),
  rootTime_       (-1.) {

  base -> options () -> GetIntegerValue ("convexification_points", nSamples_, "couenne.");

  std::string s;

  base -> options () -> GetStringValue ("convexification_type", s, "couenne.");
  if      (s == "current-point-only") convtype_ = CURRENT_ONLY;
  else if (s == "uniform-grid")       convtype_ = UNIFORM_GRID;
  else                                convtype_ = AROUND_CURPOINT;

  base -> options() -> GetStringValue("violated_cuts_only", s, "couenne."); 
  addviolated_ = (s=="yes");

  problem_ = new CouenneProblem (asl, base, jnlst_);
}


/// destructor
CouenneCutGenerator::~CouenneCutGenerator ()
  {if (problem_) delete problem_;}


/// copy constructor
CouenneCutGenerator::CouenneCutGenerator (const CouenneCutGenerator &src):

  OaDecompositionBase (src),

  firstcall_   (src. firstcall_),
  addviolated_ (src. addviolated_), 
  convtype_    (src. convtype_), 
  nSamples_    (src. nSamples_),
  problem_     (src. problem_ -> clone ()),
  nrootcuts_   (src. nrootcuts_),
  ntotalcuts_  (src. ntotalcuts_),
  septime_     (src. septime_),
  objValue_    (src. objValue_),
  nlp_         (src. nlp_),
  BabPtr_      (src. BabPtr_),
  infeasNode_  (src. infeasNode_),
  jnlst_       (src. jnlst_),
  rootTime_    (src. rootTime_) {}


#define MAX_SLOPE 1e3

/// add half-space through two points (x1,y1) and (x2,y2)
int CouenneCutGenerator::addSegment (OsiCuts &cs, int wi, int xi, 
				     CouNumber x1, CouNumber y1, 
				     CouNumber x2, CouNumber y2, int sign) const { 

  if (fabs (x2-x1) < COUENNE_EPS) {
    if (fabs (y2-y1) > MAX_SLOPE * COUENNE_EPS)
      jnlst_->Printf(J_WARNING, J_CONVEXIFYING,
		     "warning, discontinuity of %e over an interval of %e\n", y2-y1, x2-x1);
    else return createCut (cs, y2, (int) 0, wi, 1.);
  }

  CouNumber dx = x2-x1, dy = y2-y1;

  //  return createCut (cs, y1 + oppslope * x1, sign, wi, 1., xi, oppslope);
  return createCut (cs, y1*dx - dy*x1, (dx>0) ? sign : -sign, wi, dx, xi, -dy);
}


/// add tangent at (x,w) with given slope
int CouenneCutGenerator::addTangent (OsiCuts &cs, int wi, int xi, 
				     CouNumber x, CouNumber w, 
				     CouNumber slope, int sign) const
  {return createCut (cs, w - slope * x, sign, wi, 1., xi, - slope);}


/// total number of variables (original + auxiliary) of the problem
const int CouenneCutGenerator::getnvars () const
  {return problem_ -> nVars ();} 


/// Add list of options to be read from file
void CouenneCutGenerator::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> SetRegisteringCategory ("Couenne options", Bonmin::RegisteredOptions::CouenneCategory);

  roptions -> AddLowerBoundedIntegerOption
    ("convexification_cuts",
     "Specify the frequency (in terms of nodes) at which couenne ecp cuts are generated.",
     0,1,
     "A frequency of 0 amounts to never solve the NLP relaxation.");

  roptions -> AddStringOption2
    ("local_optimization_heuristic",
     "Do we search for local solutions of NLP's",
     "yes",
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
     "Apply convexification cuts before branching",
     "no",
     "no","",
     "yes","");

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_local_optimization_per_level",
     "Specify the logarithm of the number of local optimizations to perform" 
     " on average for each level of given depth of the tree.",
     -1,-1,"Solve as many nlp's at the nodes for each level of the tree. "
     "Nodes are randomly selected. If for a"
     "given level there are less nodes than this number nlp are solved for every nodes."
     "For example if parameter is 8, nlp's are solved for all node until level 8," 
     "then for half the node at level 9, 1/4 at level 10...."
     "Value -1 specify to perform at all nodes.");

  roptions -> AddStringOption3
    ("convexification_type",
     "Deterimnes in which point the linear over/under-estimator are generated",
     "current-point-only",
     "current-point-only","Only at current optimum of relaxation",
     "uniform-grid","Points chosen in a unform grid between the bounds of the problem",
     "around-current-point","At points around current optimum of relaxation");
    
  roptions -> AddLowerBoundedIntegerOption
    ("convexification_points",
     "Specify the number of points at which to convexify when convexification type"
     "is uniform-grid or arround-current-point.",
     0,1,
     "");

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

  roptions -> AddStringOption2 
    ("violated_cuts_only",
     "Yes if only violated convexification cuts should be added",
     "yes",
     "no","",
     "yes","");

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

  roptions -> setOptionExtraInfo ("branch_pt_select", 15); // Why 15? TODO

  CouenneProblem::registerOptions (roptions);
}
