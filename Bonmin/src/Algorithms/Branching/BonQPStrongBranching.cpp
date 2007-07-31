// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
// Authors: Andreas Waechter, Pierre Bonami

#include "BonminConfig.h"

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#define Verbose
#include "BonQPStrongBranching.hpp"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#include "BonBqpdSolver.hpp"
#endif

namespace Bonmin {

BonQPStrongBranching::BonQPStrongBranching(OsiTMINLPInterface * solver,
					   bool solve_nlp /* = false */) :
  BonChooseVariable(solver),
  solve_nlp_(solve_nlp)
{
}

BonQPStrongBranching::BonQPStrongBranching(const BonQPStrongBranching & rhs) :
  BonChooseVariable(rhs)
{
  solve_nlp_ = rhs.solve_nlp_;
}

BonQPStrongBranching &
BonQPStrongBranching::operator=(const BonQPStrongBranching & rhs)
{
  if (this != &rhs) {
    solve_nlp_ = rhs.solve_nlp_;
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

// Clone
ChooseVariable *
BonQPStrongBranching::clone2() const
{
  return new BonQPStrongBranching(*this);
}

int
BonQPStrongBranching::fill_changes(OsiSolverInterface * solver,
				   OsiBranchingInformation *info,
				   bool fixVariables, int numStrong,
				   double* change_down,
				   double* change_up, int& best_way)
{
  // Create a QP (or NLP) problem based on the current solution
  OsiTMINLPInterface* tminlp_interface =
    dynamic_cast<OsiTMINLPInterface*> (solver);
  TMINLP2TNLP* tminlp2tnlp = tminlp_interface->problem();

  SmartPtr<TMINLP2TNLP> branching_tqp;
  if (solve_nlp_) {
    branching_tqp = new TMINLP2TNLP(*tminlp2tnlp);
  }
  else {
    branching_tqp = new BranchingTQP(*tminlp2tnlp);
  }

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

  const Number large_number = COIN_DBL_MAX;

#ifdef OLDOLD
  bool first_solve = true;
  SmartPtr<TNLPSolver> tqp_solver =
    tminlp_interface->solver()->clone();
  // Get a warm start object from the node just solved
  CoinWarmStart* warmStart = NULL;
#ifdef COIN_HAS_FILTERSQP
  FilterSolver* filter_solver =
    dynamic_cast<FilterSolver*> (tminlp_interface->solver());
  if (filter_solver) {
    warmStart = filter_solver->getWarmStart(tminlp2tnlp);
  }
#endif
#endif //OLDOLD
  
  bool first_solve = true;
  CoinWarmStart* warmStart = NULL;
  SmartPtr<TNLPSolver> tqp_solver;
#ifdef COIN_HAS_FILTERSQP
  FilterSolver* filter_solver =
    dynamic_cast<FilterSolver*> (tminlp_interface->solver());
  if (filter_solver) {
    if (solve_nlp_) {
      warmStart = filter_solver->getWarmStart(tminlp2tnlp);
    }else {
      tqp_solver = new BqpdSolver(tminlp_interface->solver()->RegOptions(),
				  tminlp_interface->solver()->Options(),
				  tminlp_interface->solver()->Jnlst());
    }
  }
#endif
  if (IsNull(tqp_solver)) {
    tqp_solver = tminlp_interface->solver()->clone();
  }
  
  for (int i=0; i<numStrong; i++) {
    int& index = list_[i];
    const OsiObject * object = solver->object(index);
    int col_number = object->columnNumber();
    DBG_ASSERT(col_number != -1);
    const double& cutoff = info->cutoff_;

    // up
    Number curr_bnd = branching_tqp->x_l()[col_number];
    const Number up_bnd = Min(b_U[col_number],ceil(solution[col_number]));

    branching_tqp->SetVariableLowerBound(col_number, up_bnd);

    // Solve Problem
    tqp_solver->enableWarmStart();
    if (warmStart) {
      tqp_solver->setWarmStart(warmStart, branching_tqp);
    }
    TNLPSolver::ReturnStatus retstatus;
    if (first_solve) {
      retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
    }
    else {
      retstatus = tqp_solver->ReOptimizeTNLP(GetRawPtr(branching_tqp));
    }

    if (retstatus == TNLPSolver::solvedOptimal ||
	retstatus == TNLPSolver::solvedOptimalTol) {
      const Number& new_obj = branching_tqp->obj_value();
      if (new_obj > cutoff) {
	if (bb_log_level_>3) {
	  printf("Declaring branch infeasible because new_obj = %e and cutoff = %e\n",new_obj, cutoff);
	}
	delete warmStart;
	best_way = 1;
	return i;
      }
      change_up[i] = new_obj - curr_obj;
      first_solve = false;
    }
    else if (retstatus == TNLPSolver::provenInfeasible) {
      delete warmStart;
      best_way = 1;
      return i;
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
    if (warmStart) {
      tqp_solver->setWarmStart(warmStart, branching_tqp);
    }
    if (first_solve) {
      retstatus = tqp_solver->OptimizeTNLP(GetRawPtr(branching_tqp));
    }
    else {
      retstatus = tqp_solver->ReOptimizeTNLP(GetRawPtr(branching_tqp));
    }

    if (retstatus == TNLPSolver::solvedOptimal ||
	retstatus == TNLPSolver::solvedOptimalTol) {
      const Number& new_obj = branching_tqp->obj_value();
      if (new_obj > cutoff) {
	if (bb_log_level_>3) {
	  printf("Declaring branch infeasible because new_obj = %e and cutoff = %e\n",new_obj, cutoff);
	}
	delete warmStart;
	best_way = 0;
	return i;
      }
      change_down[i] = new_obj - curr_obj;
      first_solve = false;
    }
    else if (retstatus == TNLPSolver::provenInfeasible) {
      delete warmStart;
      best_way = 0;
      return i;
    }
    else {
      change_down[i] = -large_number;
    }

    branching_tqp->SetVariableUpperBound(col_number, curr_bnd);
  }

  delete warmStart;

  return -1; // nothing infeasible detected
}

}
