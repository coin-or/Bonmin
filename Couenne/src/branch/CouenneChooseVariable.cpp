/*
 * Name:    CouenneChooseVariable.cpp
 * Author:  Pietro Belotti
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneChooseVariable.hpp>

/// Default Constructor 
CouenneChooseVariable::CouenneChooseVariable () {

}

/// Constructor from solver (so we can set up arrays etc)
CouenneChooseVariable::CouenneChooseVariable (const OsiSolverInterface *solver) {

}

/// Copy constructor 
CouenneChooseVariable::CouenneChooseVariable (const CouenneChooseVariable &source):
  OsiChooseVariable (source) {

}

/// Assignment operator 
CouenneChooseVariable & CouenneChooseVariable::operator= (const CouenneChooseVariable& rhs)
{}

/// Clone
OsiChooseVariable *CouenneChooseVariable::clone() const
  {return new CouenneChooseVariable (*this);}

/// Destructor 
CouenneChooseVariable::~CouenneChooseVariable () {

}

/** Sets up strong list and clears all if initialize is true.
    Returns number of infeasibilities. 
    If returns -1 then has worked out node is infeasible!
*/

int CouenneChooseVariable::setupList (OsiBranchingInformation *info, bool initialize) {
  return 0;
}


/** Choose a variable
    Returns:
    -1 Node is infeasible
    0  Normal termination - we have a candidate
    1  All looks satisfied - no candidate
    2  We can change the bound on a variable - but we also have a strong branching candidate
    3  We can change the bound on a variable - but we have a non-strong branching candidate
    4  We can change the bound on a variable - no other candidates
    We can pick up branch from bestObjectIndex() and bestWhichWay()
    We can pick up a forced branch (can change bound) from firstForcedObjectIndex() 
    and firstForcedWhichWay()
    If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
    If fixVariables is true then 2,3,4 are all really same as problem changed
*/

int CouenneChooseVariable::chooseVariable (OsiSolverInterface * solver, 
					   OsiBranchingInformation *info, 
					   bool fixVariables) {
  return 0;
}


/// Given a candidate fill in useful information e.g. estimates
void CouenneChooseVariable::updateInformation (const OsiBranchingInformation *info,
					       int branch, OsiHotInfo * hotInfo) {

}

























#if 0



#ifdef UseOurOwn
// Initialize
// PB: ToDo Check
int 
BonChooseVariable::setupList ( OsiBranchingInformation *info, bool initialize)
{
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
  int bestPriority=INT_MAX;
  // pretend one strong even if none
  //int maximumStrong = numberStrong_ ? CoinMin(numberStrong_,numberObjects) : 1;
  int maximumStrong = numberObjects;
  int putOther = numberObjects;
  int numViolatedAtBestPriority = 0;
  int i;

  OsiObject ** object = info->solver_->objects();
  for ( i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      int priorityLevel = object[i]->priority();
      // Better priority? Flush choices.
      if (priorityLevel<bestPriority) {
	  for (int j = numViolatedAtBestPriority - 1; j >= 0; --j) {
	      list_[--putOther] = list_[j];
	      useful_[putOther] = useful_[j]; 
	  }
	  bestPriority = priorityLevel;
	  numViolatedAtBestPriority = 0;
	  check=0.0;
      } 
      if (priorityLevel==bestPriority && value>check) {
	  //add to list
	  if (numViolatedAtBestPriority < maximumStrong) {
	      list_[numViolatedAtBestPriority] = i;
	      useful_[numViolatedAtBestPriority] = value;
	      ++numViolatedAtBestPriority;
	  } else {
	      assert (useful_[checkIndex] == check);
	      list_[--putOther] = list_[checkIndex];
	      useful_[putOther] = check; 
	      list_[checkIndex] = i;
	      useful_[checkIndex] = value;
	  }
	  if (numViolatedAtBestPriority == maximumStrong) {
	      // find worst
	      check=useful_[0];
              checkIndex = 0;
	      for (int j = 1; j < maximumStrong; ++j) {
		  if (useful_[j] < check) {
		      check = useful_[j];
		      checkIndex = j;
		  }
	      }
	  }
      } else {
	  // to end
	  list_[--putOther] = i;
	  useful_[putOther] = value;
      }
    }
  }
  // Get list
  numberOnList_ = numViolatedAtBestPriority;
  if (numberOnList_) {
      for (i = 0; i < numberOnList_; ++i) {
	  useful_[i] = - useful_[i];
      }
      // Sort 
      CoinSort_2(useful_,useful_+numberOnList_,list_);
      // move others
      i = numberOnList_;
      for (;putOther<numberObjects;putOther++, i++) {
	  list_[i]=list_[putOther];
	  useful_[i] = - useful_[putOther];
      }
      assert (i==numberUnsatisfied_);
      if (!numberStrong_)
	  numberOnList_=0;
  }

  if (bb_log_level_>4) {
    for (int i=0; i<Min(numberUnsatisfied_,numberStrong_); i++)
      printf("list_[%5d] = %5d, usefull_[%5d] = %23.16e\n", i,list_[i],i,useful_[i]);
  }
  return numberUnsatisfied_;
}
#endif
 
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
int 
BonChooseVariable::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
  if (numberUnsatisfied_) {

    int numStrong = Min(numberUnsatisfied_, numberStrong_);

    const Number large_number = COIN_DBL_MAX;
 
    int best_i = 0;
    int best_way;
    double best_change = large_number;
    if (numStrong > 1) {
      // space for storing the predicted changes in the objective
      double* change_up = new double[numStrong];
      double* change_down = new double[numStrong];

      // This will return -1 unless an infeasible node was detected
      best_i = fill_changes(solver, info, fixVariables,
			    numStrong, change_down, change_up, best_way);
      if (best_i == -1) {
	// Determine most promising branching variable
	best_i = -1;
	best_change = -large_number;
	for (int i=0; i<numStrong; i++) {
	  if (bb_log_level_>4) {
	    printf("i = %d down = %15.6e up = %15.6e\n",
		   i, change_down[i], change_up[i]);
	  }

	  // for now, we look for the best combined change
	  double change_min = Min(change_down[i], change_up[i]);
	  double change_max = Max(change_down[i], change_up[i]);
	  double change_comp = 2.*change_max + change_min;
	  // only use new value if significantly larger (rel_fact)
	  const Number rel_fact = 1e-6;
	  if (best_change - change_comp <
	      -rel_fact*Max(1.,fabs(best_change))) {
	    best_change = change_comp;
	    best_i = i;
	  }
	}

	if (best_i == -1) {
	  best_i = 0;
	}
	assert(best_i != -1);
      }
      delete [] change_up;
      delete [] change_down;
    }

    if (bb_log_level_>3) {
      printf("best_i = %d  best_change = %15.6e\n", best_i, best_change);
    }

    bestObjectIndex_=list_[best_i];
    bestWhichWay_ = solver->object(bestObjectIndex_)->whichWay();
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    return 0;
  } else {
    return 1;
  }
}

}
#endif
