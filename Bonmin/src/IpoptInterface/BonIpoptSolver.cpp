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

namespace Bonmin{
  IpoptSolver::~IpoptSolver(){}

  ///virtual constructor
  TNLPSolver * 
  IpoptSolver::createNew()
  { return new IpoptSolver();}


  void
  IpoptSolver::Initialize(std::string params_file)
  {
    app_.Initialize(params_file);
    setMinlpDefaults(app_.Options());
  }

  void
  IpoptSolver::Initialize(std::istream &is)
  {
    app_.Initialize(is);
    setMinlpDefaults(app_.Options());
  }

  TNLPSolver::ReturnStatus
  IpoptSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
    TNLPSolver::ReturnStatus ret_status;
    if(!zeroDimension(tnlp, ret_status))
      {
	optimizationStatus_ = app_.OptimizeTNLP(tnlp);
      }
    else
      {
	if(ret_status == solvedOptimal)
	  optimizationStatus_ = Ipopt::Solve_Succeeded;
	else if(ret_status == provenInfeasible)
	  optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
      }
    return solverReturnStatus(optimizationStatus_);
  }


  TNLPSolver::ReturnStatus
  IpoptSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp)
  {
    TNLPSolver::ReturnStatus ret_status;
    if(!zeroDimension(tnlp, ret_status))
      {
	optimizationStatus_ = app_.ReOptimizeTNLP(tnlp);
      }
    else
      {
	if(ret_status == solvedOptimal)
	  optimizationStatus_ = Ipopt::Solve_Succeeded;
	else if(ret_status == provenInfeasible)
	  optimizationStatus_ = Ipopt::Infeasible_Problem_Detected;
      }
    return solverReturnStatus(optimizationStatus_);
  }

  Ipopt::SmartPtr<Ipopt::RegisteredOptions>
  IpoptSolver::RegOptions()
  {
    return app_.RegOptions();
  }

  Ipopt::SmartPtr<const Ipopt::OptionsList>
  IpoptSolver::Options() const
  {
    return app_.Options();
  }

  Ipopt::SmartPtr<Ipopt::OptionsList>
  IpoptSolver::Options()
  {
    return app_.Options();
  }

  /// Get the CpuTime of the last optimization.
  double 
  IpoptSolver::CPUTime()
  {
    const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_.Statistics();
    if(IsValid(stats))
      {
	return stats->TotalCPUTime();
      }
    else
      {
	throw CoinError("No statistics available from Ipopt","CPUTime","Bonmin::IpoptSolver");
      }
  }

  /// Get the iteration count of the last optimization.
  int
  IpoptSolver::IterationCount()
  {
    const Ipopt::SmartPtr<Ipopt::SolveStatistics>  stats = app_.Statistics();
    if(IsValid(stats))
      {
	return stats->IterationCount();
      }
    else
      {
	throw CoinError("No statistics available from Ipopt","IterationCount","Bonmin::IpoptSolver");
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
      return computationError;
    case Ipopt::Not_Enough_Degrees_Of_Freedom:
    case Ipopt::Invalid_Problem_Definition:
      return illDefinedProblem;
    case Ipopt::Invalid_Option:
    case Ipopt::Invalid_Number_Detected:
      return illegalOption;
    case Ipopt::NonIpopt_Exception_Thrown:
      return externalException;
    case Ipopt::Insufficient_Memory:
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
}
