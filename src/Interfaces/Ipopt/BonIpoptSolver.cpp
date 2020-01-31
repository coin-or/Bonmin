// (C) Copyright International Business Machines (IBM) 2005, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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

namespace Bonmin
{

  std::string IpoptSolver::solverName_ = "Ipopt";

  /// Constructor
  IpoptSolver::IpoptSolver(bool createEmpty /*= false*/):
      TNLPSolver(),
      problemHadZeroDimension_(false),
      warmStartStrategy_(1),
      enable_warm_start_(false),
      optimized_before_(false)
  {
    if (createEmpty) return;
    app_ = new Ipopt::IpoptApplication(GetRawPtr(roptions_), options_, journalist_);
#ifdef NO_CATCH_ALL
    app_->RethrowNonIpoptException(true);
#endif
  }

  /// Constructor with Passed in journalist, registered options, options
  IpoptSolver::IpoptSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
      Ipopt::SmartPtr<Ipopt::OptionsList> options,
      Ipopt::SmartPtr<Ipopt::Journalist> journalist,
      const std::string & prefix):
      TNLPSolver(roptions, options, journalist, prefix),
      problemHadZeroDimension_(false),
      warmStartStrategy_(1),
      enable_warm_start_(false),
      optimized_before_(false)
  {
    roptions_ = roptions;
    app_ = new Ipopt::IpoptApplication(GetRawPtr(roptions), options, journalist);
#ifdef NO_CATCH_ALL
    app_->RethrowNonIpoptException(true);
#endif
  }

  /// Constructor with Passed in journalist, registered options, options
  IpoptSolver::IpoptSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
      Ipopt::SmartPtr<Ipopt::OptionsList> options,
      Ipopt::SmartPtr<Ipopt::Journalist> journalist):
      TNLPSolver(roptions, options, journalist, "bonmin."),
      problemHadZeroDimension_(false),
      warmStartStrategy_(1),
      enable_warm_start_(false),
      optimized_before_(false)
  {
    roptions_ = roptions;
    app_ = new Ipopt::IpoptApplication(GetRawPtr(roptions), options, journalist);
#ifdef NO_CATCH_ALL
    app_->RethrowNonIpoptException(true);
#endif
  }

  IpoptSolver::~IpoptSolver()
  {}

  IpoptSolver::IpoptSolver(const IpoptSolver &other):
    TNLPSolver(other),
    optimizationStatus_(other.optimizationStatus_),
    problemHadZeroDimension_(other.problemHadZeroDimension_),
    warmStartStrategy_(other.warmStartStrategy_),
    enable_warm_start_(false),
    optimized_before_(false){
      app_ = new Ipopt::IpoptApplication(GetRawPtr(roptions_), options_, journalist_);
#ifdef NO_CATCH_ALL
      app_->RethrowNonIpoptException(true);
#endif
  }

  ///virtual constructor
  Ipopt::SmartPtr<TNLPSolver>
  IpoptSolver::clone()
  {
    Ipopt::SmartPtr<IpoptSolver> retval = new IpoptSolver(*this);
    retval->app_->Initialize("");
    retval->default_log_level_ = default_log_level_;
    return GetRawPtr(retval);
  }


  bool
  IpoptSolver::Initialize(std::string params_file)
  {
    Ipopt::ApplicationReturnStatus status =
      app_->Initialize(params_file);
    if (status != Ipopt::Solve_Succeeded) {
      return false;
    }
    options_->GetEnumValue("warm_start",warmStartStrategy_,prefix());
    setMinlpDefaults(options_);
    optimized_before_ = false;
    return true;
  }

  bool
  IpoptSolver::Initialize(std::istream &is)
  {
    Ipopt::ApplicationReturnStatus status =
      app_->Initialize(is);
    if (status != Ipopt::Solve_Succeeded) {
      return false;
    }
    options_->GetEnumValue("warm_start",warmStartStrategy_,prefix());
    setMinlpDefaults(app_->Options());
    optimized_before_ = false;
    return true;
  }

  TNLPSolver::ReturnStatus
  IpoptSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
#if 0
    printf("Global Time limit set to %g\n", time_limit_);
    double local_time_limit = time_limit_ - 
                              CoinCpuTime() + start_time_;
    printf("Time limit set to %g\n", local_time_limit);
    if(local_time_limit <= 0.){
       optimizationStatus_ = Ipopt::Maximum_CpuTime_Exceeded;
       return solverReturnStatus(optimizationStatus_);
    }
#endif
    TNLPSolver::ReturnStatus ret_status;
    if (!zeroDimension(tnlp, ret_status)) {
#if 0
      if(time_limit_ < DBL_MAX){
           options_->SetNumericValue("max_cpu_time", local_time_limit,
                                  true, true);
      }
#endif
      if (enable_warm_start_ && optimized_before_) {
        optimizationStatus_ = app_->ReOptimizeTNLP(tnlp);
      }
      else {
        optimizationStatus_ = app_->OptimizeTNLP(tnlp);
      }
      optimized_before_ = true;
      problemHadZeroDimension_ = false;
    }
    else {
      problemHadZeroDimension_ = true;
      if (ret_status == solvedOptimal)
        optimizationStatus_ = Ipopt::Solve_Succeeded;
      else if (ret_status == provenInfeasible)
        optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
    }

    if (BonminAbortAll)
      optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;

    return solverReturnStatus(optimizationStatus_);
  }


  TNLPSolver::ReturnStatus
  IpoptSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
#if 0
    printf("Global Time limit set to %g\n", time_limit_);
    double local_time_limit = time_limit_ - 
                              CoinCpuTime() + start_time_;
    printf("Time limit set to %g\n", local_time_limit);
    if(local_time_limit <= 0.){
       optimizationStatus_ = Ipopt::Maximum_CpuTime_Exceeded;
       return solverReturnStatus(optimizationStatus_);
    }
#endif
    TNLPSolver::ReturnStatus ret_status;
    if (!zeroDimension(tnlp, ret_status)) {
#if 0
      if(time_limit_ < DBL_MAX){
        options_->SetNumericValue("max_cpu_time", 
                                  std::max(0., local_time_limit),
                                  true, true);
      }
#endif
      if (optimized_before_) {
        optimizationStatus_ = app_->ReOptimizeTNLP(tnlp);
      }
      else {
        optimizationStatus_ = app_->OptimizeTNLP(tnlp);
      }
      problemHadZeroDimension_ = false;
      optimized_before_ = true;
    }
    else {
      problemHadZeroDimension_ = true;
      if (ret_status == solvedOptimal)
        optimizationStatus_ = Ipopt::Solve_Succeeded;
      else if (ret_status == provenInfeasible)
        optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
    }
    if (BonminAbortAll)
      optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
    return solverReturnStatus(optimizationStatus_);
  }

  /// Get the CpuTime of the last optimization.
  double
  IpoptSolver::CPUTime()
  {
    if (problemHadZeroDimension_) {
      return 0.;
    }
    else {
      const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_->Statistics();
      if (IsValid(stats)) {
        return stats->TotalCpuTime();
      }
      else {
        app_->Jnlst()->Printf(Ipopt::J_WARNING, Ipopt::J_STATISTICS, "No statistics available from Ipopt in Bonmin::IpoptSolver::CPUTime\n");
        return 0.;
        //throw CoinError("No statistics available from Ipopt","CPUTime","Bonmin::IpoptSolver");
      }
    }
  }

  /// Get the iteration count of the last optimization.
  int
  IpoptSolver::IterationCount()
  {
    if (problemHadZeroDimension_) {
      return 0;
    }
    else {
      const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_->Statistics();
      if (IsValid(stats)) {
        return stats->IterationCount();
      }
      else {
        app_->Jnlst()->Printf(Ipopt::J_WARNING, Ipopt::J_STATISTICS, "No statistics available from Ipopt in Bonmin::IpoptSolver::IterationCount\n");
        return 0;
        //throw CoinError("No statistics available from Ipopt","IterationCount","Bonmin::IpoptSolver");
      }

    }
  }


  void
  IpoptSolver::setMinlpDefaults(Ipopt::SmartPtr<Ipopt::OptionsList> Options)
  {
    bool set = false;
    double dummy_dbl;
    int dummy_int;
    set = Options->GetNumericValue("gamma_phi", dummy_dbl, "");
    if(!set)
    Options->SetNumericValue("gamma_phi", 1e-8, true, true);
    set = Options->GetNumericValue("gamma_theta", dummy_dbl, "");
    if(!set)
    Options->SetNumericValue("gamma_theta", 1e-4, true, true);
    set = Options->GetNumericValue("required_infeasibility_reduction", dummy_dbl, "");
    if(!set)
    Options->SetNumericValue("required_infeasibility_reduction", 0.1, true, true);
    set = Options->GetEnumValue("expect_infeasible_problem",dummy_int, "");
    if(!set)
    Options->SetStringValue("expect_infeasible_problem","yes", true, true);
    set = Options->GetEnumValue("mu_strategy", dummy_int, "");
    if(!set)
    Options->SetStringValue("mu_strategy", "adaptive", true, true);
    set = Options->GetEnumValue("mu_oracle",dummy_int, "");
    if(!set)
    Options->SetStringValue("mu_oracle","probing", true, true);
    if(!Options->GetIntegerValue("print_level",default_log_level_,"")) {
      default_log_level_ = 1;
      Options->SetIntegerValue("print_level",1, true, true);
    }
  }


  ////////////////////////////////////////////////////////////////////
  // Methods returning info on how the solution process terminated  //
  ////////////////////////////////////////////////////////////////////
  /// Are there a numerical difficulties?
  TNLPSolver::ReturnStatus IpoptSolver::solverReturnStatus(Ipopt::ApplicationReturnStatus optimization_status) const
  {

    switch (optimization_status) {
    case Ipopt::Maximum_Iterations_Exceeded:
    case Ipopt::User_Requested_Stop:
    case Ipopt::Restoration_Failed:
      return iterationLimit;
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
    case Ipopt::Feasible_Point_Found:
      return solvedOptimal;
    case Ipopt::Search_Direction_Becomes_Too_Small:
      return doesNotConverge;
    case Ipopt::Solved_To_Acceptable_Level:
      return solvedOptimalTol;
    case Ipopt::Infeasible_Problem_Detected:
      return provenInfeasible;
    case Ipopt::Diverging_Iterates:
      return unbounded;
    case Ipopt::Maximum_CpuTime_Exceeded:
      return timeLimit;
    default:
      return exception;
    }
  }

/// Get warm start used in last optimization
CoinWarmStart *
IpoptSolver::getUsedWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const
{
  if(tnlp->x_init() == NULL || tnlp->duals_init() == NULL)
    return NULL;
  return  new IpoptWarmStart(tnlp->num_variables(),
                             2*tnlp->num_variables() + 
                             tnlp->num_constraints(),
                             tnlp->x_init(), tnlp->duals_init());
}
/// Get warmstarting information
  CoinWarmStart*
  IpoptSolver::getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const
  {
      if (warmStartStrategy_==2) {
        Ipopt::SmartPtr<IpoptInteriorWarmStarter> warm_starter =
          Ipopt::SmartPtr<IpoptInteriorWarmStarter>(tnlp->GetWarmStarter());
        return new IpoptWarmStart(tnlp, warm_starter);
      }
      else  return new IpoptWarmStart(tnlp, NULL);
  }


  bool
  IpoptSolver::setWarmStart(const CoinWarmStart* warmstart,
      Ipopt::SmartPtr<TMINLP2TNLP> tnlp)
  {
    if (!warmstart && warmStartStrategy_)
      return 0;
    const IpoptWarmStart * ws = dynamic_cast<const IpoptWarmStart*> (warmstart);
    if(ws == NULL) return 0;
    if (ws->empty())//reset initial point and leave
    {
      disableWarmStart();
      return 1;
    }
    if(ws->dualSize() > 0){
      tnlp->setDualsInit(ws->dualSize(), ws->dual());
      enableWarmStart();
    }
    else
      disableWarmStart();
#ifndef NDEBUG
    int numcols = tnlp->num_variables();
    int numrows = tnlp->num_constraints();
#endif

    assert(numcols == ws->primalSize());
    assert(2*numcols + numrows == ws->dualSize());
    tnlp->setxInit(ws->primalSize(), ws->primal());

    if (IsValid(ws->warm_starter()))
      tnlp->SetWarmStarter(ws->warm_starter());

    return 1;
  }

  bool 
  IpoptSolver::warmStartIsValid(const CoinWarmStart * ws) const{
    const IpoptWarmStart* iws = dynamic_cast<const IpoptWarmStart*>(ws);
    if (iws && !iws->empty()) {
      return true;
    }
    return false;
  }

  CoinWarmStart *
  IpoptSolver::getEmptyWarmStart() const
  {
    return new IpoptWarmStart(1);
  }

  void
  IpoptSolver::enableWarmStart()
  {
    enable_warm_start_ = true;
    options_->SetStringValue("warm_start_init_point", "yes");
  }

  void
  IpoptSolver::disableWarmStart()
  {
    enable_warm_start_ = false;
    options_->SetStringValue("warm_start_init_point", "no");
  }


  void
  IpoptSolver::setOutputToDefault()
  {
     options_->SetIntegerValue("print_level", default_log_level_, true, true);
  }


  void
  IpoptSolver::forceSolverOutput(int log_level)
  {
     options_->SetIntegerValue("print_level", log_level, true, true);
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
    if (errorNum() >=0)
      return errorNames[errorNum()];
    if (errorNum() == -1) return errorNames[6];
    else if (errorNum() == -2) return errorNames[7];
    else if (errorNum() == -3) return errorNames[8];
    else if (errorNum() == -10) return errorNames[9];
    else if (errorNum() == -11) return errorNames[10];
    else if (errorNum() == -12) return errorNames[11];
    else if (errorNum() == -13) return errorNames[12];
    else if (errorNum() == -100) return errorNames[13];
    else if (errorNum() == -101) return errorNames[14];
    else if (errorNum() == -102) return errorNames[15];
    else if (errorNum() == -199) return errorNames[16];
    throw CoinError("UnsolvedError", "UnsolvedError::errorName()","Unrecognized optimization status in ipopt.");
  }

  std::string IpoptSolver::UnsolvedIpoptError::solverName_ = "Ipopt";

  const std::string &
  IpoptSolver::UnsolvedIpoptError::solverName() const
  {
    return solverName_;
  }




}//end namespace Bonmin
