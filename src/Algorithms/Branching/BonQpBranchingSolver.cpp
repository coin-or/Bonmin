// Copyright (C) 2006, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.
// Authors: Andreas Waechter, Pierre Bonami

#include "BonminConfig.h"

#include "CoinPragma.hpp"

#define Verbose
#include "BonQpBranchingSolver.hpp"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#include "BonBqpdSolver.hpp"
#endif

namespace Bonmin
{

  QpBranchingSolver::QpBranchingSolver(OsiTMINLPInterface * solver)
      :
      StrongBranchingSolver(solver)
  {}

  QpBranchingSolver::QpBranchingSolver(const QpBranchingSolver & rhs) :
      StrongBranchingSolver(rhs)
  {}

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
#ifdef TIME_BQPD
     printf("QPBRANCH Timings for %i sbs\n", times_.numsolve);
     printf("QPBRANCH %i pivots\n", times_.pivots);
     printf("QPBRANCH Creating  : %g\n", times_.create);
     printf("QPBRANCH Solving : %g\n", times_.solve);
     printf("QPBRANCH Warming  : %g\n", times_.warm_start);
     printf("QPBRANCH Resolving : %g\n", times_.resolve);
#endif
  }

  void QpBranchingSolver::
  markHotStart(OsiTMINLPInterface* tminlp_interface)
  {
    TMINLP2TNLP* tminlp2tnlp = tminlp_interface->problem();
    branching_tqp_ = new BranchingTQP(tminlp2tnlp);

    first_solve_ = true;
#ifdef COIN_HAS_FILTERSQP
    FilterSolver* filter_solver =
      dynamic_cast<FilterSolver*> (tminlp_interface->solver());
    if (filter_solver) {
      Ipopt::SmartPtr<BqpdSolver> qp_solver =
	new BqpdSolver(RegOptions(), Options(), Jnlst(), tminlp_interface->prefix());
#if 1
      // Solve the QP with the original bounds and set the hot start
      // information
      TNLPSolver::ReturnStatus retstatus;
      retstatus = qp_solver->OptimizeTNLP(GetRawPtr(branching_tqp_));
      if (retstatus == TNLPSolver::solvedOptimal ||
	  retstatus == TNLPSolver::solvedOptimalTol) {
	first_solve_ = false;
	qp_solver->markHotStart();
      }
#endif
      tqp_solver_ = GetRawPtr(qp_solver);
      //tqp_solver_ = new FilterSolver(RegOptions(), Options(), Jnlst());
    }
#endif
    if (IsNull(tqp_solver_)) {
      tqp_solver_ = tminlp_interface->solver()->clone();
    }
    tqp_solver_->enableWarmStart();
  }

  void QpBranchingSolver::
  unmarkHotStart(OsiTMINLPInterface* tminlp_interface)
  {
    // Free memory
#ifdef TIME_BQPD
    BqpdSolver * qp_solver = dynamic_cast<BqpdSolver *>(GetRawPtr(tqp_solver_));
    if(qp_solver) times_ += qp_solver->times();
#endif
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
