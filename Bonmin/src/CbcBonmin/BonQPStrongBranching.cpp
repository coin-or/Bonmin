// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonQPStrongBranching.hpp"

namespace Bonmin {

BonQPStrongBranching::BonQPStrongBranching(OsiTMINLPInterface * solver) :
  BonChooseVariable(solver)
{
}

BonQPStrongBranching::BonQPStrongBranching(const BonQPStrongBranching & rhs) :
  BonChooseVariable(rhs)
{
}

BonQPStrongBranching &
BonQPStrongBranching::operator=(const BonQPStrongBranching & rhs)
{
  if (this != &rhs) {
    BonChooseVariable::operator=(rhs);
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
#define RunAllProblems
int 
BonQPStrongBranching::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
#ifdef RunAllProblems
  if (fixVariables) {
    chooseVariable(solver,info, false);

    printf("== Beginning: Strong branching with curvature estimator:\n");
    BonChooseVariable::chooseVariable(solver,info, fixVariables);

    printf("== Beginning: Strong branching with QP:\n");
  }
  else {
    printf("== Beginning: Strong branching with NLP:\n");
  }
#endif

  if (numberUnsatisfied_) {

    // Create a QP problem based on the current solution
    OsiTMINLPInterface* tminlp_interface =
      dynamic_cast<OsiTMINLPInterface*> (solver);
    const TMINLP2TNLP* tminlp2tnlp = tminlp_interface->problem();

#ifdef RunAllProblems
    SmartPtr<TMINLP2TNLP> branching_tqp;
    if (fixVariables) {
      branching_tqp = new BranchingTQP(*tminlp2tnlp);
    }
    else {
      branching_tqp = new TMINLP2TNLP(*tminlp2tnlp);
    }
#else
    //SmartPtr<BranchingTQP> branching_tqp = new BranchingTQP(*tminlp2tnlp);
    SmartPtr<TMINLP2TNLP> branching_tqp = new TMINLP2TNLP(*tminlp2tnlp);
#endif
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
    int best_way = -1;
    if (numStrong > 1) {
      bool first_solve = true;
      SmartPtr<TNLPSolver> tqp_solver =
	tminlp_interface->solver()->clone();
      //TODO: activate this: tqp_solver->enableWarmStart();

      for (int i=0; i<numStrong; i++) {
	int& index = list_[i];
	const OsiObject * object = solver->object(index);
	int col_number = object->columnNumber();
	DBG_ASSERT(col_number != -1);

	// up
	Number curr_bnd = branching_tqp->x_l()[col_number];
	const Number up_bnd = Min(b_U[col_number],ceil(solution[col_number]));

#ifdef VeryVerbose
	printf("up bounds: %d sol %e cur %e new %e\n", col_number, solution[col_number], curr_bnd, up_bnd);
#endif
	branching_tqp->SetVariableLowerBound(col_number, up_bnd);

	// Solve Problem
	TNLPSolver::ReturnStatus retstatus;
	if (first_solve) {
	  retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	}
	else {
	  retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	  //retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
	}
#ifdef VeryVerbose
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
#ifdef VeryVerbose
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
	  printf("i = %d down = %15.6e up = %15.6e\n", i,change_down[i], change_up[i]);
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
    printf("best_i = %d  best_change = %15.6e best_way = %d\n", best_i, best_change, best_way);
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
