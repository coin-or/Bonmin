// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonLpBranchingSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include <vector>

namespace Bonmin {

LpBranchingSolver::LpBranchingSolver(OsiTMINLPInterface * solver) :
  StrongBranchingSolver(solver),
  lin_(NULL),
  warm_(NULL),
  ecp_(NULL)
{
  SmartPtr<TNLPSolver> tnlp_solver =
    static_cast<TNLPSolver *> (solver->solver());
  SmartPtr<OptionsList> options = tnlp_solver->Options();

  options->GetIntegerValue("ecp_max_rounds_strong",
			   maxCuttingPlaneIterations_, "bonmin.");
  options->GetNumericValue("ecp_abs_tol_strong", abs_ecp_tol_,"bonmin.");
  options->GetNumericValue("ecp_rel_tol_strong", rel_ecp_tol_,"bonmin.");
}

LpBranchingSolver::LpBranchingSolver(const LpBranchingSolver & rhs) :
  StrongBranchingSolver(rhs),
  lin_(NULL),
  warm_(NULL),
  ecp_(NULL),
  maxCuttingPlaneIterations_(rhs.maxCuttingPlaneIterations_)
{
}

LpBranchingSolver &
LpBranchingSolver::operator=(const LpBranchingSolver & rhs)
{
  if (this != &rhs) {
    StrongBranchingSolver::operator=(rhs);
  }
  maxCuttingPlaneIterations_ = rhs.maxCuttingPlaneIterations_;
  abs_ecp_tol_ = rhs.abs_ecp_tol_;
  rel_ecp_tol_ = rhs.rel_ecp_tol_;
  // I assume that no LP solver information is ever copied
  delete lin_;
  delete warm_;
  delete ecp_;
  lin_ = NULL;
  warm_ = NULL;
  ecp_ = NULL;
  return *this;
}

LpBranchingSolver::~LpBranchingSolver ()
{
  delete lin_;
  delete warm_;
  delete ecp_; 
}

void LpBranchingSolver::
markHotStart(OsiTMINLPInterface* tminlp_interface)
{
  lin_ = new OsiClpSolverInterface();
  tminlp_interface->extractLinearRelaxation(*lin_, true, false);
  // ???lin_->setDblParam(OsiDualObjectiveLimit, info->cutoff_);
  lin_->messageHandler()->setLogLevel(0);
  lin_->resolve();
  warm_ = lin_->getWarmStart();
  if (maxCuttingPlaneIterations_)
    ecp_ = new EcpCuts(tminlp_interface, maxCuttingPlaneIterations_,
		       abs_ecp_tol_, rel_ecp_tol_, -1.);
}

void LpBranchingSolver::
unmarkHotStart(OsiTMINLPInterface* tminlp_interface)
{
  // Free memory
  delete lin_;
  delete warm_;
  delete ecp_;
  lin_ = NULL;
  warm_ = NULL;
  ecp_ = NULL;
}

TNLPSolver::ReturnStatus LpBranchingSolver::
solveFromHotStart(OsiTMINLPInterface* tminlp_interface)
{
  TNLPSolver::ReturnStatus retstatus = TNLPSolver::iterationLimit;

  // updated the bounds of the linear solver
  std::vector<int> diff_low_bnd_index;
  std::vector<double> diff_low_bnd_value;
  std::vector<int> diff_up_bnd_index;
  std::vector<double> diff_up_bnd_value;

#if 0
  // deleteme
  for (int i=0; i<tminlp_interface->getNumCols(); i++) {
    printf("%3d ol = %e nl = %e   ou = %e nu = %e\n",i,tminlp_interface->getColLower()[i],lin_->getColLower()[i],tminlp_interface->getColUpper()[i],lin_->getColUpper()[i]);
  }
#endif
  // Get the bounds.  We assume that the bounds in the linear solver
  // are always the original ones
  const int numCols = tminlp_interface->getNumCols();
  const double* colLow_orig = lin_->getColLower();
  const double* colUp_orig = lin_->getColUpper();
  const double* colLow = tminlp_interface->getColLower();
  const double* colUp = tminlp_interface->getColUpper();
  // Set the bounds on the LP solver according to the changes in
  // tminlp_interface
  for (int i=0; i<numCols; i++) {
    const double& lo = colLow[i];
    if (colLow_orig[i] < lo) {
      diff_low_bnd_value.push_back(colLow_orig[i]);
      diff_low_bnd_index.push_back(i);
      lin_->setColLower(i,lo);
    }
    const double& up = colUp[i];
    if (colUp_orig[i] > up) {
      diff_up_bnd_index.push_back(i);
      diff_up_bnd_value.push_back(colUp_orig[i]);
      lin_->setColUpper(i,lo);
    }
  }

#if 0
  // deleteme
  for (int i=0; i<numCols; i++) {
    printf("%3d ol = %e nl = %e   ou = %e nu = %e\n",i,tminlp_interface->getColLower()[i],lin_->getColLower()[i],tminlp_interface->getColUpper()[i],lin_->getColUpper()[i]);
  }
#endif

  lin_->setWarmStart(warm_);
  lin_->resolve();

  double obj = lin_->getObjValue();
  bool go_on = true;
  if (lin_->isProvenPrimalInfeasible()) {
    retstatus = TNLPSolver::provenInfeasible;
    go_on = false;
  }
  else if (lin_->isIterationLimitReached()) {
    retstatus = TNLPSolver::iterationLimit;
    go_on = false;
  }
  else {
    if (maxCuttingPlaneIterations_ > 0 && go_on) {
      double violation;
      obj = ecp_->doEcpRounds(*lin_, true, &violation);
      if (obj == COIN_DBL_MAX) {
	retstatus = TNLPSolver::provenInfeasible;
      }
      else if (violation <= 1e-8) {
	retstatus = TNLPSolver::solvedOptimal;
      }
    }
  }
  tminlp_interface->problem()->set_obj_value(obj);

  //restore the original bounds
  for (unsigned int i = 0; i < diff_low_bnd_index.size(); i++) {
    lin_->setColLower(diff_low_bnd_index[i],diff_low_bnd_value[i]);
  }
  for (unsigned int i = 0; i < diff_up_bnd_index.size(); i++) {
    lin_->setColUpper(diff_up_bnd_index[i],diff_up_bnd_value[i]);
  }

  return retstatus;
}

void
LpBranchingSolver::registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
  roptions->AddLowerBoundedIntegerOption
    ("ecp_max_rounds_strong",
     "Set the maximal number of rounds of ECP cuts in strong branching.",
     0,5,
     "");
  roptions->AddLowerBoundedNumberOption
    ("ecp_abs_tol_strong",
     "Set the absolute termination tolerance for ECP rounds in strong branching.",
     0,false,1e-6,
     "");
  roptions->AddLowerBoundedNumberOption
    ("ecp_rel_tol_strong",
     "Set the relative termination tolerance for ECP rounds in strong branching.",
     0,false,1e-1,
     "");
}

}
