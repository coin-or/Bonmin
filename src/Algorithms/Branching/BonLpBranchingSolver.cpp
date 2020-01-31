// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "BonLpBranchingSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include <vector>

namespace Bonmin
{

  LpBranchingSolver::LpBranchingSolver(BabSetupBase * b) :
      StrongBranchingSolver(b->nonlinearSolver()),
      lin_(NULL),
      warm_(NULL),
      ecp_(NULL)
  {
    Ipopt::SmartPtr<TNLPSolver> tnlp_solver =
       static_cast<TNLPSolver *> (b->nonlinearSolver()->solver());
    Ipopt::SmartPtr<Ipopt::OptionsList> options = tnlp_solver->options();

	    options->GetIntegerValue("ecp_max_rounds_strong",
	                             maxCuttingPlaneIterations_,
                                     b->nonlinearSolver()->prefix());
	    options->GetNumericValue("ecp_abs_tol_strong",
                                     abs_ecp_tol_,
                                     b->nonlinearSolver()->prefix());
	    options->GetNumericValue("ecp_rel_tol_strong",
                                     rel_ecp_tol_,
                                     b->nonlinearSolver()->prefix());
	    int dummy;
	    options->GetEnumValue("lp_strong_warmstart_method", 
                                  dummy,
                                  b->nonlinearSolver()->prefix());
	    warm_start_mode_ = (WarmStartMethod) dummy;
	  }

  LpBranchingSolver::LpBranchingSolver(const LpBranchingSolver & rhs) :
      StrongBranchingSolver(rhs),
      lin_(NULL),
      warm_(NULL),
      ecp_(NULL),
      maxCuttingPlaneIterations_(rhs.maxCuttingPlaneIterations_),
      abs_ecp_tol_(rhs.abs_ecp_tol_),
      rel_ecp_tol_(rhs.rel_ecp_tol_),
      warm_start_mode_(rhs.warm_start_mode_)
  {}

  LpBranchingSolver &
  LpBranchingSolver::operator=(const LpBranchingSolver & rhs)
  {
    if (this != &rhs) {
      StrongBranchingSolver::operator=(rhs);
    }
    maxCuttingPlaneIterations_ = rhs.maxCuttingPlaneIterations_;
    abs_ecp_tol_ = rhs.abs_ecp_tol_;
    rel_ecp_tol_ = rhs.rel_ecp_tol_;
    warm_start_mode_ = rhs.warm_start_mode_;
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
    tminlp_interface->extractLinearRelaxation(*lin_, tminlp_interface->getColSolution(),
                                              true);
    double cutoff = -DBL_MAX;
    tminlp_interface->getDblParam(OsiDualObjectiveLimit, cutoff);
    lin_->setDblParam(OsiDualObjectiveLimit, cutoff);
    //printf("Cutoff %g # ecp iteration %i\n",cutoff, maxCuttingPlaneIterations_);
    lin_->messageHandler()->setLogLevel(0);
    lin_->resolve();
    warm_ = lin_->getWarmStart();
    //if (maxCuttingPlaneIterations_)
    //  ecp_ = new EcpCuts(tminlp_interface, maxCuttingPlaneIterations_,
    //      abs_ecp_tol_, rel_ecp_tol_, -1.);
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
    TNLPSolver::ReturnStatus retstatus = TNLPSolver::solvedOptimal;

    // updated the bounds of the linear solver
    std::vector<int> diff_low_bnd_index;
    std::vector<double> diff_low_bnd_value;
    std::vector<int> diff_up_bnd_index;
    std::vector<double> diff_up_bnd_value;

    // Get the bounds.  We assume that the bounds in the linear solver
    // are always the original ones
    const int numCols = tminlp_interface->getNumCols();
    const double* colLow_orig = lin_->getColLower();
    const double* colUp_orig = lin_->getColUpper();
    const double* colLow = tminlp_interface->getColLower();
    const double* colUp = tminlp_interface->getColUpper();

    OsiSolverInterface * lin = lin_;
    // eventualy clone lin_
    if(warm_start_mode_ == Clone){
      lin = lin_->clone();
//      std::cout<<"Cloning it"<<std::endl;
    }
    // Set the bounds on the LP solver according to the changes in
    // tminlp_interface
    for (int i=0; i<numCols; i++) {
      const double& lo = colLow[i];
      if (colLow_orig[i] < lo) {
        if(warm_start_mode_ == Basis){
          diff_low_bnd_value.push_back(colLow_orig[i]);
          diff_low_bnd_index.push_back(i);
        }
        lin->setColLower(i,lo);
      }
      const double& up = colUp[i];
      if (colUp_orig[i] > up) {
        if(warm_start_mode_ == Basis){
          diff_up_bnd_index.push_back(i);
          diff_up_bnd_value.push_back(colUp_orig[i]);
        }
        lin->setColUpper(i,lo);
      }
    }

    if(warm_start_mode_ == Basis){
      lin->setWarmStart(warm_);
    }

    lin->resolve();

    double obj = lin->getObjValue();
    bool go_on = true;
    if (lin->isProvenPrimalInfeasible() || 
        lin->isDualObjectiveLimitReached()) {
      retstatus = TNLPSolver::provenInfeasible;
      go_on = false;
    }
    else if (lin->isIterationLimitReached()) {
      retstatus = TNLPSolver::iterationLimit;
      go_on = false;
    }
    else {
      if (maxCuttingPlaneIterations_ > 0 && go_on) {
        double violation;
        obj = ecp_->doEcpRounds(*lin, true, &violation);
        if (obj == COIN_DBL_MAX) {
          retstatus = TNLPSolver::provenInfeasible;
        }
        else if (violation <= 1e-8) {
          retstatus = TNLPSolver::solvedOptimal;
        }
      }
    }
    tminlp_interface->problem()->set_obj_value(obj);
    tminlp_interface->problem()->Set_x_sol(numCols, lin_->getColSolution());

    //restore the original bounds
    if(warm_start_mode_ == Basis){
      for (unsigned int i = 0; i < diff_low_bnd_index.size(); i++) {
        lin_->setColLower(diff_low_bnd_index[i],diff_low_bnd_value[i]);
      }
      for (unsigned int i = 0; i < diff_up_bnd_index.size(); i++) {
        lin_->setColUpper(diff_up_bnd_index[i],diff_up_bnd_value[i]);
      }
    }
    else {
      delete lin;
    }
    return retstatus;
  }

  void
  LpBranchingSolver::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("ECP based strong branching",RegisteredOptions::UndocumentedCategory);
    roptions->AddLowerBoundedIntegerOption
    ("ecp_max_rounds_strong",
     "Set the maximal number of rounds of ECP cuts in strong branching.",
     0,0,
     "");
    roptions->setOptionExtraInfo("ecp_max_rounds_strong",63);
    roptions->AddLowerBoundedNumberOption
    ("ecp_abs_tol_strong",
     "Set the absolute termination tolerance for ECP rounds in strong branching.",
     0,false,1e-6,
     "");
    roptions->setOptionExtraInfo("ecp_abs_tol_strong",63);
    roptions->AddLowerBoundedNumberOption
    ("ecp_rel_tol_strong",
     "Set the relative termination tolerance for ECP rounds in strong branching.",
     0,false,1e-1,
     "");
    roptions->setOptionExtraInfo("ecp_rel_tol_strong",63);
    roptions->AddStringOption2
    ("lp_strong_warmstart_method",
     "Choose method to use for warm starting lp in strong branching",
     "Basis",
     "Basis", "Use optimal basis of node",
     "Clone", "Clone optimal problem of node",
     "(Advanced stuff)");
    roptions->setOptionExtraInfo("lp_strong_warmstart_method",63);
  }

}
