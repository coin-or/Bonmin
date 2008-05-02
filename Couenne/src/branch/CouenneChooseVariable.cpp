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

/// Clone
OsiChooseVariable *CouenneChooseVariable::clone() const
{return new CouenneChooseVariable (*this);}

/// Destructor 
CouenneChooseVariable::~CouenneChooseVariable () {}


#ifndef CCS_EXPERIMENTAL

/** Sets up strong list and clears all if initialize is true.
    Returns number of infeasibilities. 
    If returns -1 then has worked out node is infeasible!
*/

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
#endif

// Returns true if solution looks feasible against given objects
bool CouenneChooseVariable::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects)

{return problem_ -> checkNLP (solution, solution [problem_ -> Obj (0) -> Body () -> Index ()]);}



/// Add list of options to be read from file
void CouenneChooseVariable::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> AddStringOption2
    ("branch_fbbt",
     "Apply bound tightening before branching",
     "yes",
     "no","",
     "yes","");

  roptions -> AddStringOption2
    ("branch_conv_cuts",
     "Apply convexification cuts before branching",
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
}


#ifdef CCS_EXPERIMENTAL
// Initialize
int CouenneChooseVariable::setupList (OsiBranchingInformation *info, 
				      bool initialize) {

  if (initialize) {
    status_=-2;
    delete [] goodSolution_;
    bestObjectIndex_=-1;
    numberStrongDone_=0;
    numberStrongIterations_ = 0;
    numberStrongFixed_ = 0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
  }
  numberOnList_=0;
  numberUnsatisfied_=0;
  int numberObjects = solver_->numberObjects();
  assert (numberObjects);
  double check = 0.0;
  int checkIndex=0;
  int bestPriority=COIN_INT_MAX;
  // pretend one strong even if none
  int maximumStrong= numberStrong_ ? CoinMin(numberStrong_,numberObjects) : 1;
  int putOther = numberObjects;
  int i;
  for (i=0;i<maximumStrong;i++) {
    list_[i]=-1;
    useful_[i]=0.0;
  }
  OsiObject ** object = info->solver_->objects();
  // Say feasible
  bool feasible = true;
  //printf ("setup loop %d objs-----\n", numberObjects);
  for ( i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      if (value==COIN_DBL_MAX) {
	// infeasible
	feasible=false;
	break;
      }
      int priorityLevel = object[i]->priority();
      // Better priority? Flush choices.
      if (priorityLevel<bestPriority) {
	for (int j=0;j<maximumStrong;j++) {
	  if (list_[j]>=0) {
	    int iObject = list_[j];
	    list_[j]=-1;
	    useful_[j]=0.0;
	    list_[--putOther]=iObject;
	  }
	}
	bestPriority = priorityLevel;
	check=0.0;
      } 
      if (priorityLevel==bestPriority) {
	if (value>check) {
	  //add to list
	  int iObject = list_[checkIndex];
	  if (iObject>=0)
	    list_[--putOther]=iObject;  // to end
	  list_[checkIndex]=i;
	  useful_[checkIndex]=value;
	  // find worst
	  check=COIN_DBL_MAX;
	  for (int j=0;j<maximumStrong;j++) {
	    if (list_[j]>=0) {
	      if (useful_[j]<check) {
		check=useful_[j];
		checkIndex=j;
	      }
	    } else {
	      check=0.0;
	      checkIndex = j;
	      break;
	    }
	  }
	} else {
	  // to end
	  list_[--putOther]=i;
	}
      } else {
	// to end
	list_[--putOther]=i;
      }
    }
  }
  //printf ("--------- end setup loop\n");
  // Get list
  numberOnList_=0;
  if (feasible) {
    for (i=0;i<maximumStrong;i++) {
      if (list_[i]>=0) {
	list_[numberOnList_]=list_[i];
	useful_[numberOnList_++]=-useful_[i];
      }
    }
    if (numberOnList_) {
      // Sort 
      CoinSort_2(useful_,useful_+numberOnList_,list_);
      // move others
      i = numberOnList_;
      for (;putOther<numberObjects;putOther++) 
	list_[i++]=list_[putOther];
      assert (i==numberUnsatisfied_);
      if (!numberStrong_)
	numberOnList_=0;
    } 
  } else {
    // not feasible
    numberUnsatisfied_=-1;
  }
  return numberUnsatisfied_;
}


/* Choose a variable
   Returns - 
   -1 Node is infeasible
   0  Normal termination - we have a candidate
   1  All looks satisfied - no candidate
   2  We can change the bound on a variable - but we also have a strong branching candidate
   3  We can change the bound on a variable - but we have a non-strong branching candidate
   4  We can change the bound on a variable - no other candidates
   We can pick up branch from whichObject() and whichWay()
   We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
   If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
*/
int CouenneChooseVariable::chooseVariable (OsiSolverInterface * solver, 
					   OsiBranchingInformation *info, 
					   bool fixVariables) {
  if (numberUnsatisfied_) {
    bestObjectIndex_=list_[0];
    bestWhichWay_ = solver->object(bestObjectIndex_)->whichWay();
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    return 0;
  } else {
    return 1;
  }
}

#endif
