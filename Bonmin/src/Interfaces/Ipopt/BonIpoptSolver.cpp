// (C) Copyright International Business Machines (IBM) 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, IBM
//
// Date : 26/09/2006


#include "BonIpoptSolver.hpp"
#include "IpSolveStatistics.hpp"
#include "CoinError.hpp"

#include "BonIpoptInteriorWarmStarter.hpp"
#include "BonIpoptWarmStart.hpp"


extern bool BonminAbortAll;

namespace Bonmin{

  std::string IpoptSolver::solverName_ = "Ipopt";

  /// Constructor
  IpoptSolver::IpoptSolver(bool createEmpty /*= false*/):
    problemHadZeroDimension_(false),
    warmStartStrategy_(1),
    enable_warm_start_(false),
    optimized_before_(false)
  {
    if(createEmpty) return;
    app_ = new Ipopt::IpoptApplication();
  }

  /// Constructor with Passed in journalist, registered options, options
  IpoptSolver::IpoptSolver(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions,
                           Ipopt::SmartPtr<Ipopt::OptionsList> options,
                           Ipopt::SmartPtr<Ipopt::Journalist> journalist):
  TNLPSolver(roptions, options, journalist),
  problemHadZeroDimension_(false),
  warmStartStrategy_(1),
  enable_warm_start_(false),
  optimized_before_(false)
{
    app_ = new Ipopt::IpoptApplication(roptions, options, journalist);
}
  
  IpoptSolver::~IpoptSolver(){}

  ///virtual constructor
  Ipopt::SmartPtr<TNLPSolver> 
  IpoptSolver::clone()
  {
    SmartPtr<IpoptSolver> retval = new IpoptSolver(true);

    retval->warmStartStrategy_ = warmStartStrategy_;
    retval->problemHadZeroDimension_ = problemHadZeroDimension_;
    retval->optimizationStatus_ = optimizationStatus_;

    retval->app_ = app_->clone();
    enable_warm_start_ = false;
    optimized_before_ = false;

    return GetRawPtr(retval);
  }


  void
  IpoptSolver::Initialize(std::string params_file)
  {
    app_->Initialize(params_file);
    app_->Options()->GetEnumValue("warm_start",warmStartStrategy_,"bonmin.");
    setMinlpDefaults(app_->Options());
    optimized_before_ = false;
  }

  void
  IpoptSolver::Initialize(std::istream &is)
  {
    app_->Initialize(is);
    app_->Options()->GetEnumValue("warm_start",warmStartStrategy_,"bonmin.");
    setMinlpDefaults(app_->Options());
    optimized_before_ = false;
  }

  TNLPSolver::ReturnStatus
  IpoptSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
    TNLPSolver::ReturnStatus ret_status;
    if(!zeroDimension(tnlp, ret_status))
      {
	if (enable_warm_start_ && optimized_before_) {
	  optimizationStatus_ = app_->ReOptimizeTNLP(tnlp);
	}
	else {
	  optimizationStatus_ = app_->OptimizeTNLP(tnlp);
	}
	optimized_before_ = true;
	problemHadZeroDimension_ = false;
      }
    else
      {
	problemHadZeroDimension_ = true;
	if(ret_status == solvedOptimal)
	  optimizationStatus_ = Ipopt::Solve_Succeeded;
	else if(ret_status == provenInfeasible)
	  optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
      }

     if(BonminAbortAll)
       optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;

    return solverReturnStatus(optimizationStatus_);
  }


  TNLPSolver::ReturnStatus
  IpoptSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
    TNLPSolver::ReturnStatus ret_status;
    if(!zeroDimension(tnlp, ret_status))
      {
	if (optimized_before_) {	
	  optimizationStatus_ = app_->ReOptimizeTNLP(tnlp);
	}
	else {
	  optimizationStatus_ = app_->OptimizeTNLP(tnlp);
	}
	problemHadZeroDimension_ = false;
	optimized_before_ = true;
      }
    else
      {
	problemHadZeroDimension_ = true;
	if(ret_status == solvedOptimal)
	  optimizationStatus_ = Ipopt::Solve_Succeeded;
	else if(ret_status == provenInfeasible)
	  optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
      }
     if(BonminAbortAll)
       optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
      return solverReturnStatus(optimizationStatus_);
  }

  Ipopt::SmartPtr<Ipopt::Journalist>
  IpoptSolver::Jnlst()
  {return app_->Jnlst();}

  Ipopt::SmartPtr<Ipopt::RegisteredOptions>
  IpoptSolver::RegOptions()
  {
    return app_->RegOptions();
  }

  Ipopt::SmartPtr<const Ipopt::OptionsList>
  IpoptSolver::Options() const
  {
    return ConstPtr(app_)->Options();
  }

  Ipopt::SmartPtr<Ipopt::OptionsList>
  IpoptSolver::Options()
  {
    return app_->Options();
  }

  /// Get the CpuTime of the last optimization.
  double 
  IpoptSolver::CPUTime()
  {
    if(problemHadZeroDimension_)
      {
	return 0.;
      }
    else
      {
	const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_->Statistics();
	if(IsValid(stats))
	  {
	    return stats->TotalCPUTime();
	  }
	else
	  {
            printf("TODO: No statistics available from Ipopt in Bonmin::IpoptSolver::CPUTime\n");
            return 0.;
	    //throw CoinError("No statistics available from Ipopt","CPUTime","Bonmin::IpoptSolver");
	  }
      }
  }

  /// Get the iteration count of the last optimization.
  int
  IpoptSolver::IterationCount()
  {
    if(problemHadZeroDimension_)
      {
	return 0;
      }
    else
      {
	const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_->Statistics();
	if(IsValid(stats))
	  {
	    return stats->IterationCount();
	  }
	else
	  {
            printf("TODO: No statistics available from Ipopt in Bonmin::IpoptSolver::IterationCount\n");
            return 0;
	    //throw CoinError("No statistics available from Ipopt","IterationCount","Bonmin::IpoptSolver");
	  }
	
      }
  }


  void 
  IpoptSolver::setMinlpDefaults(Ipopt::SmartPtr<Ipopt::OptionsList> Options){
    Options->SetNumericValue("gamma_phi", 1e-8, true, true);
    Options->SetNumericValue("gamma_theta", 1e-4, true, true);
    Options->SetNumericValue("required_infeasibility_reduction", 0.1, true, true);
    Options->SetStringValue("expect_infeasible_problem","yes", true, true);
    Options->SetStringValue("mu_strategy", "adaptive", true, true);
    Options->SetStringValue("mu_oracle","probing", true, true);
    Options->SetIntegerValue("print_level",1, true, true);
  }


  ////////////////////////////////////////////////////////////////////
  // Methods returning info on how the solution process terminated  //
  ////////////////////////////////////////////////////////////////////
  /// Are there a numerical difficulties?
  TNLPSolver::ReturnStatus IpoptSolver::solverReturnStatus(Ipopt::ApplicationReturnStatus optimization_status) const
  {
  
    switch(optimization_status){
    case Ipopt::Maximum_Iterations_Exceeded:
    case Ipopt::User_Requested_Stop:
      return iterationLimit;
    case Ipopt::Restoration_Failed:
    case Ipopt::Error_In_Step_Computation:
    case Ipopt::Unrecoverable_Exception:
    case Ipopt::Insufficient_Memory:
      return computationError;
    case Ipopt::Not_Enough_Degrees_Of_Freedom:
      return notEnoughFreedom;
    case Ipopt::Invalid_Problem_Definition:
      return illDefinedProblem;
    case Ipopt::Invalid_Option:
    case Ipopt::Invalid_Number_Detected:
      return illegalOption;
    case Ipopt::NonIpopt_Exception_Thrown:
      return externalException;
    case Ipopt::Internal_Error:
      return exception;
    case Ipopt::Solve_Succeeded:
      return solvedOptimal;
    case Ipopt::Search_Direction_Becomes_Too_Small:
      std::cerr<<"Warning : need to verify that Search_Direction_Becomes_Too_Small is indeed Ok"<<std::endl;
    case Ipopt::Solved_To_Acceptable_Level:
      return solvedOptimalTol;
    case Ipopt::Infeasible_Problem_Detected:
      return provenInfeasible;
    case Ipopt::Diverging_Iterates:
      return unbounded;
    default:
      return exception;
    }
  }

/// Get warmstarting information
CoinWarmStart*
IpoptSolver::getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const
{
  if(warmStartStrategy_) {
    if(warmStartStrategy_==2) {
        Ipopt::SmartPtr<IpoptInteriorWarmStarter> warm_starter =
        Ipopt::SmartPtr<IpoptInteriorWarmStarter>(tnlp->GetWarmStarter());
        return new IpoptWarmStart(tnlp, warm_starter);
    }
    else  return new IpoptWarmStart(tnlp, NULL);
  }
  else
     return new IpoptWarmStart(tnlp->num_variables(), tnlp->num_constraints());
}


bool
IpoptSolver::setWarmStart(const CoinWarmStart* warmstart,
                          Ipopt::SmartPtr<TMINLP2TNLP> tnlp)
{
  if(!warmstart && warmStartStrategy_)
    return 0;
  const IpoptWarmStart * ws = dynamic_cast<const IpoptWarmStart*> (warmstart);
  if(ws->empty())//reset initial point and leave
  {
    disableWarmStart();
    return 1;
  }
  enableWarmStart();
  int numcols = tnlp->num_variables();
  int numrows = tnlp->num_constraints();
  const double * colLow = tnlp->x_l();
  const double * colUp = tnlp->x_u();
  for(int i = 0; i < numcols ; i++) {
    CoinWarmStartBasis::Status status = ws->getStructStatus(i);
    if(status == CoinWarmStartBasis::atLowerBound) {
      tnlp->setxInit(i,colLow[i]);
      tnlp->setDualInit(i + numcols + numrows,0.);
    }
    else if(status == CoinWarmStartBasis::atUpperBound) {
      tnlp->setxInit(i,colUp[i]);
      tnlp->setDualInit(i + numrows,0.);
    }
    else {
      tnlp->setDualInit(i + numrows,0.);
      tnlp->setDualInit(i + numcols + numrows, 0.);
    }
  }
  for(int i = 0; i < numrows ; i++) {
    CoinWarmStartBasis::Status status = ws->getArtifStatus(i);
    if(status == CoinWarmStartBasis::atLowerBound) {
      tnlp->setDualInit(i,0.);
    }
  }
  int nElem = ws->values()->getNumElements();
  const int * inds = ws->values()->getIndices();
  const double * elems = ws->values()->getElements();

  for(int i = 0 ; i < nElem ; i++) {
    tnlp->setxInit(inds[i],elems[i]);
  }

  if(IsValid(ws->warm_starter()))
    tnlp->SetWarmStarter(ws->warm_starter());
  return 1;
}


CoinWarmStart * 
IpoptSolver::getEmptyWarmStart() const
{return new IpoptWarmStart(1);}

  void 
  IpoptSolver::enableWarmStart()
  {
    enable_warm_start_ = true;
    Options()->SetStringValue("warm_start_init_point", "yes");
  }

  void 
  IpoptSolver::disableWarmStart()
  {
    enable_warm_start_ = false;
    Options()->SetStringValue("warm_start_init_point", "no");
  }


  void
  IpoptSolver::turnOffOutput(){
  }


  void
  IpoptSolver::turnOnOutput(){
  }


/*******************************************************************************/
// Class for throwing errors reported from Ipopt
/******************************************************************************/

std::string
IpoptSolver::UnsolvedIpoptError::errorNames[17] ={"Solve succeeded",
    "Solved to acceptable level",
    "Infeasible problem detected",
    "Search direction becomes too small",
    "Diverging iterates",
    "User requested stop",
    "Maximum iterations exceeded",
    "Restoration failed",
    "Error in step computation",
    "Not enough degrees of freedom",
    "Invalid problem definition",
    "Invalid option",
    "Invalid number detected",
    "Unrecoverable exception",
    "NonIpopt exception thrown",
    "Insufficient memory",
    "Internal error"};

const std::string &
IpoptSolver::UnsolvedIpoptError::errorName() const
{
  if(errorNum() >=0)
    return errorNames[errorNum()];
  if(errorNum() == -1) return errorNames[6];
  else if(errorNum() == -2) return errorNames[7];
  else if(errorNum() == -3) return errorNames[8];
  else if(errorNum() == -10) return errorNames[9];
  else if(errorNum() == -11) return errorNames[10];
  else if(errorNum() == -12) return errorNames[11];
  else if(errorNum() == -13) return errorNames[12];
  else if(errorNum() == -100) return errorNames[13];
  else if(errorNum() == -101) return errorNames[14];
  else if(errorNum() == -102) return errorNames[15];
  else if(errorNum() == -199) return errorNames[16];
  throw CoinError("UnsolvedError", "UnsolvedError::errorName()","Unrecognized optimization status in ipopt.");
}

std::string IpoptSolver::UnsolvedIpoptError::solverName_ = "Ipopt";

const std::string &
IpoptSolver::UnsolvedIpoptError::solverName() const
{
   return solverName_;
}




}//end namespace Bonmin
