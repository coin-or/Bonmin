// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
// Authors: Andreas Waechter, Pierre Bonami

#include "BonminConfig.h"

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#define Verbose
#include "BonQpBranchingSolver.hpp"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#include "BonBqpdSolver.hpp"
#endif

namespace Bonmin {

QpBranchingSolver::QpBranchingSolver(OsiTMINLPInterface * solver)
  :
  StrongBranchingSolver(solver)
{
}

QpBranchingSolver::QpBranchingSolver(const QpBranchingSolver & rhs) :
  StrongBranchingSolver(rhs)
{
}

QpBranchingSolver &
QpBranchingSolver::operator=(const QpBranchingSolver & rhs)
{
  if (this != &rhs) {
    StrongBranchingSolver::operator=(rhs);
  }
  return *this;
}

QpBranchingSolver::~QpBranchingSolver ()
{
}

void QpBranchingSolver::
markHotStart(OsiTMINLPInterface* tminlp_interface)
{
  TMINLP2TNLP* tminlp2tnlp = tminlp_interface->problem();
  branching_tqp_ = new BranchingTQP(tminlp2tnlp);

#ifdef COIN_HAS_FILTERSQP
  FilterSolver* filter_solver =
    dynamic_cast<FilterSolver*> (tminlp_interface->solver());
  if (filter_solver) {
    //tqp_solver_ = new BqpdSolver(RegOptions(), Options(), Jnlst());
    tqp_solver_ = new FilterSolver(RegOptions(), Options(), Jnlst());
  }
#endif
  if (IsNull(tqp_solver_)) {
    tqp_solver_ = tminlp_interface->solver()->clone();
  }
  tqp_solver_->enableWarmStart();
  first_solve_ = true;
}

void QpBranchingSolver::
unmarkHotStart(OsiTMINLPInterface* tminlp_interface)
{
  // Free memory
  branching_tqp_ = NULL;
  tqp_solver_ = NULL;
}

TNLPSolver::ReturnStatus QpBranchingSolver::
solveFromHotStart(OsiTMINLPInterface* tminlp_interface)
{
  TNLPSolver::ReturnStatus retstatus;
  if (first_solve_) {
    retstatus = tqp_solver_->OptimizeTNLP(GetRawPtr(branching_tqp_));
  }
  else {
    retstatus = tqp_solver_->ReOptimizeTNLP(GetRawPtr(branching_tqp_));
  }

  if (retstatus == TNLPSolver::solvedOptimal ||
      retstatus == TNLPSolver::solvedOptimalTol) {
    // don't way we solve the problem, since otherwise the pseudo costs
    // are updated and that is maybe not so good???
    //retstatus = TNLPSolver::iterationLimit;
    first_solve_ = false;
  }
  //retstatus = TNLPSolver::iterationLimit;
  return retstatus;
}

}
