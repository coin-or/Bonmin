// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonChooseVariable.hpp"
#include "IpBlas.hpp"

namespace Bonmin {

BonChooseVariable::BonChooseVariable(OsiTMINLPInterface * solver) :
  OsiChooseVariable(solver)
{
  SmartPtr<TNLPSolver> tnlp_solver =
    dynamic_cast<TNLPSolver *> (solver->solver());
  SmartPtr<Journalist> jnlst = tnlp_solver->Jnlst();
  SmartPtr<OptionsList> options = tnlp_solver->Options();
  SmartPtr<TNLP> tnlp = solver->problem();

  cur_estimator_ = new CurvatureEstimator(jnlst, options, tnlp);
}

BonChooseVariable::BonChooseVariable(const BonChooseVariable & rhs) :
  OsiChooseVariable(rhs)
{
  cur_estimator_ = rhs.cur_estimator_;
}

BonChooseVariable &
BonChooseVariable::operator=(const BonChooseVariable & rhs)
{
  if (this != &rhs) {
    OsiChooseVariable::operator=(rhs);
    cur_estimator_ = rhs.cur_estimator_;
  }
  return *this;
}

BonChooseVariable::~BonChooseVariable ()
{
}

// Clone
OsiChooseVariable *
BonChooseVariable::clone() const
{
  return new BonChooseVariable(*this);
}

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
  int maximumStrong= Min(numberStrong_ ? numberStrong_ : 1, numberObjects);
  int putOther = numberObjects;
  int i;
  for (i=0;i<maximumStrong;i++) {
    list_[i]=-1;
    useful_[i]=0.0;
  }
  OsiObject ** object = info->solver_->objects();
  for ( i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
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
      }
    }
  }
  // Get list
  numberOnList_=0;
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
int 
BonChooseVariable::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
  if (numberUnsatisfied_) {

    // Get info about the current solution
    const double * solution = solver->getColSolution();// Current solution
    int numCols = solver->getNumCols();
    int numRows = solver->getNumRows();
    const double * lam = solver->getRowPrice();
    const double * z_L = lam + numRows;
    const double * z_U = z_L + numCols;
    const double * b_L = solver->getColLower();
    const double * b_U = solver->getColUpper();

    int numberObjects = solver_->numberObjects();
    int numStrong = Min(numberObjects, numberStrong_);
    // space for storing the predicted changes in the objective
    double * change_up = new double[numStrong];
    double * change_down = new double[numStrong];

    // Set up stuff for the curvature estimator
    bool new_bounds = true;
    bool new_x = true;
    bool new_mults = true;
    double * orig_d = new double[numCols];
    double * projected_d = new double[numCols];
    const Number zero = 0.;
    IpBlasDcopy(numCols, &zero, 0, orig_d, 1);

    const Number large_number = COIN_DBL_MAX;
 
    for (int i=0; i<numStrong; i++) {
      int index = list_[i];
      const OsiObject * object = solver->object(index);
      int col_number = object->columnNumber();
      DBG_ASSERT(col_number != -1);
      // up
      orig_d[col_number] = 1.;
      double gradLagTd;
      double dTHLagd;
      bool retval =
	cur_estimator_->ComputeNullSpaceCurvature(
          new_bounds, numCols, solution, new_x, z_L, z_U,
	  numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);
      // ToDo
      assert(retval);
      new_bounds = false;
      new_x = false;
      new_mults = false;

      // Determine step size and predicted change
      const double &curr_val = solution[col_number];
      if (retval) {
	const double up_val = Min(b_U[col_number],ceil(curr_val));
	double alpha = (up_val-curr_val)/projected_d[col_number];
	change_up[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
      }
      else {
	change_up[i] = large_number;
      }

      // down
      orig_d[col_number] = -1.;
      retval = cur_estimator_->ComputeNullSpaceCurvature(
          new_bounds, numCols, solution, new_x, z_L, z_U,
	  numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);

      // Determine step size and predicted change
      if (retval) {
	const double down_val = Max(b_L[col_number],floor(curr_val));
	double alpha = (down_val-curr_val)/projected_d[col_number];
	change_down[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
      }
      else {
	change_down[i] = large_number;
      }

      orig_d[col_number] = 0.;
    }
    delete [] orig_d;
    delete [] projected_d;

    // Determine most promising branching variable
    int best_i = -1;
    double best_change = large_number;
    for (int i=0; i<numStrong; i++) {
      // for now, we look for the best combined change
      double change_min = Min(change_down[i], change_up[i]);
      double change_max = Max(change_down[i], change_up[i]);
      double change_comp = 2.*change_max + change_min;
      if (best_change > change_comp) {
	best_change = change_comp;
	best_i = i;
      }
    }

    delete [] change_up;
    delete [] change_down;

    assert(best_i != -1);

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
