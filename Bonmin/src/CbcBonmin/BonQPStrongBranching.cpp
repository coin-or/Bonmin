// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonQPStrongBranching.hpp"

namespace Bonmin {

BonQPStrongBranching::BonQPStrongBranching(OsiTMINLPInterface * solver) :
  OsiChooseVariable(solver)
{
}

BonQPStrongBranching::BonQPStrongBranching(const BonQPStrongBranching & rhs) :
  OsiChooseVariable(rhs)
{
}

BonQPStrongBranching &
BonQPStrongBranching::operator=(const BonQPStrongBranching & rhs)
{
  if (this != &rhs) {
    OsiChooseVariable::operator=(rhs);
  }
  return *this;
}

BonQPStrongBranching::~BonQPStrongBranching ()
{
}

// Clone
OsiChooseVariable *
BonQPStrongBranching::clone() const
{
  return new BonQPStrongBranching(*this);
}

#define Verbose
// For now there is no difference to what John has, so let's just use his
#ifdef UseOurOwn
// Initialize
// PB: ToDo Check
int 
BonQPStrongBranching::setupList ( OsiBranchingInformation *info, bool initialize)
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
  int maximumStrong = numberStrong_ ? CoinMin(numberStrong_,numberObjects) : 1;
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
  //  DELETEME
#ifdef Verbose
  printf("numberOnList_: %i, numberUnsatisfied_: %i, numberStrong_: %i \n",
	 numberOnList_, numberUnsatisfied_, numberStrong_);
  for (int i=0; i<Min(numberUnsatisfied_,numberStrong_); i++)
    printf("list_[%5d] = %5d, usefull_[%5d] = %23.16e\n", i,list_[i],i,useful_[i]);
#endif
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
BonQPStrongBranching::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
  if (numberUnsatisfied_) {

    // Create a QP problem based on the current solution
    OsiTMINLPInterface* tminlp_interface =
      dynamic_cast<OsiTMINLPInterface*> (solver);
    const TMINLP2TNLP* tminlp2tnlp = tminlp_interface->problem();
    SmartPtr<BranchingTQP> branching_tqp = new BranchingTQP(*tminlp2tnlp);
    const Number curr_obj = tminlp2tnlp->obj_value();

    // Get info about the current solution
    const double * solution = solver->getColSolution();// Current solution
    //    int numCols = solver->getNumCols();
    //    int numRows = solver->getNumRows();
    //    const double * lam = solver->getRowPrice();
    //    const double * z_L = lam + numRows;
    //    const double * z_U = z_L + numCols;
    const double * b_L = solver->getColLower();
    const double * b_U = solver->getColUpper();

    int numStrong = Min(numberUnsatisfied_, numberStrong_);
    // space for storing the predicted changes in the objective
    double * change_up = new double[numStrong];
    double * change_down = new double[numStrong];

    const Number large_number = COIN_DBL_MAX;
 
    int best_i = 0;
    double best_change = large_number;
    int found_infeasible = -1;
    int best_way;
    if (numStrong > 1) {
      bool first_solve = true;
      SmartPtr<TNLPSolver> tqp_solver =
	tminlp_interface->solver()->clone();

      for (int i=0; i<numStrong; i++) {
	int& index = list_[i];
	const OsiObject * object = solver->object(index);
	int col_number = object->columnNumber();
	DBG_ASSERT(col_number != -1);

	// up
	Number curr_bnd = branching_tqp->x_l()[col_number];
	const Number up_bnd = Min(b_U[col_number],ceil(solution[col_number]));

#ifdef Verbose
	printf("up bounds: %d sol %e cur %e new %e\n", col_number, solution[col_number], curr_bnd, up_bnd);
#endif
	branching_tqp->SetVariableLowerBound(col_number, up_bnd);

	// Solve Problem
	TNLPSolver::ReturnStatus retstatus;
	if (first_solve) {
	  retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	}
	else {
	  retstatus = tqp_solver->ReOptimizeTNLP(GetRawPtr(branching_tqp));
	  //retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	}
#ifdef Verbose
	// DELETEME
	printf("up: retstatus = %d obj = %e\n", retstatus, branching_tqp->obj_value());
#endif
	
	if (retstatus == TNLPSolver::solvedOptimal ||
	    retstatus == TNLPSolver::solvedOptimalTol) {
	  change_up[i] = branching_tqp->obj_value() - curr_obj;
	  first_solve = false;
	}
	else if (retstatus == TNLPSolver::provenInfeasible) {
	  // We try this for now - we should probably skip the rest of
	  // the tests
	  found_infeasible = i;
	  best_way = 1;
	  break;
	  change_up[i] = 0.;//large_number;
	  first_solve = false;
	}
	else {
	  change_up[i] = -large_number;
	}

	branching_tqp->SetVariableLowerBound(col_number, curr_bnd);

	// down
	curr_bnd = branching_tqp->x_u()[col_number];
	const Number down_bnd = Max(b_L[col_number],floor(solution[col_number]));
	branching_tqp->SetVariableUpperBound(col_number, down_bnd);

	// Solve Problem
	if (first_solve) {
	  retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	}
	else {
	  retstatus = tqp_solver->ReOptimizeTNLP(GetRawPtr(branching_tqp));
	}
#ifdef Verbose
	// DELETEME
	printf("down: retstatus = %d obj = %e\n", retstatus, branching_tqp->obj_value());
#endif

	if (retstatus == TNLPSolver::solvedOptimal ||
	    retstatus == TNLPSolver::solvedOptimalTol) {
	  change_down[i] = branching_tqp->obj_value() - curr_obj;
	  first_solve = false;
	}
	else if (retstatus == TNLPSolver::provenInfeasible) {
	  found_infeasible = i;
	  best_way = 0;
	  break;
	  // We try this for now - we should probably skip the rest of
	  // the tests
	  change_down[i] = 0.;//large_number;
	  first_solve = false;
	}
	else {
	  change_down[i] = -large_number;
	}

	branching_tqp->SetVariableUpperBound(col_number, curr_bnd);
      }

      // Determine most promising branching variable
      best_i = -1;
      best_change = -large_number;
      if (found_infeasible>-1) {
	best_i = found_infeasible;
	best_change = large_number;
      }
      else {
	for (int i=0; i<numStrong; i++) {
#ifdef Verbose
	//DELETEME
	  printf("i = %d down = %e up = %e\n", i,change_down[i], change_up[i]);
#endif
	  // for now, we look for the best combined change
	  double change_min = Min(change_down[i], change_up[i]);
	  double change_max = Max(change_down[i], change_up[i]);
	  double change_comp = 2.*change_max + change_min;
	  // only use new value if significantly larger (rel_fact)
	  const Number rel_fact = 1e-6;
	  if (best_change*(1.+rel_fact) < change_comp) {
	    best_change = change_comp;
	    best_i = i;
	    if (change_down[i] < change_up[i]) {
	      best_way = 1;
	    }
	    else {
	      best_way = 0;
	    }
	  }
	}
      }

      if (best_i == -1) {
	best_i = 0;
      }
      assert(best_i != -1);
    }
    delete [] change_up;
    delete [] change_down;

#ifdef Verbose
    //DELETEME
    printf("best_i = %d  best_change = %e best_way = %d\n", best_i, best_change, best_way);
#endif

    bestObjectIndex_=list_[best_i];
    bestWhichWay_ = best_way; // AW: check! solver->object(bestObjectIndex_)->whichWay();
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    return 0;
  } else {
    return 1;
  }
}

}
