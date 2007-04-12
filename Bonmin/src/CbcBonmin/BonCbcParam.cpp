// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#include "BonOsiTMINLPInterface.hpp"
#include "BonCbcParam.hpp"

namespace Bonmin
{
  bool
  BonminCbcParam::extractParams(OsiTMINLPInterface * solver)
  {
    bool success = true;

    Ipopt::SmartPtr<Ipopt::OptionsList> Options = solver->retrieve_options();

    //extract IpoptInterface special params
    solver->extractInterfaceParams();

    //log levels
    Options->GetIntegerValue("bb_log_level",bbLogLevel,"bonmin.");
    Options->GetIntegerValue("bb_log_interval",logInterval,"bonmin.");
    Options->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    Options->GetIntegerValue("milp_log_level",milpLogLevel,"bonmin.");
    Options->GetIntegerValue("oa_log_level",oaLogLevel,"bonmin.");
    Options->GetNumericValue("oa_log_frequency",oaLogFrequency,"bonmin.");
    Options->GetIntegerValue("nlp_log_level",nlpLogLevel,"bonmin.");
    //General options
    Options->GetEnumValue("algorithm",algo,"bonmin.");
    Options->GetNumericValue("time_limit", maxTime, "bonmin.");
    Options->GetIntegerValue("node_limit",maxNodes,"bonmin.");
    Options->GetIntegerValue("solution_limit",maxSolutions,"bonmin.");
    Options->GetIntegerValue("iteration_limit",maxIterations,"bonmin.");
    Options->GetNumericValue("integer_tolerance",intTol,"bonmin.");
    Options->GetNumericValue("allowable_gap",allowableGap,"bonmin.");
    Options->GetNumericValue("allowable_fraction_gap",allowableFractionGap,"bonmin.");
    Options->GetNumericValue("cutoff_decr",cutoffDecr,"bonmin.");
    Options->GetNumericValue("cutoff",cutoff,"bonmin.");

    // Branch & bound setting
    Options->GetEnumValue("nodeselect_stra",nodeSelection,"bonmin.");
    Options->GetEnumValue("varselect_stra",varSelection,"bonmin.");
    Options->GetIntegerValue("number_strong_branch",numberStrong,"bonmin.");
    Options->GetIntegerValue("number_before_trust", minReliability,"bonmin.");
    Options->GetIntegerValue("number_ecp_rounds_strong", numEcpRoundsStrong,"bonmin.");

    Options->GetEnumValue("sos_constraints",disableSos,"bonmin.");
    // Robustness and non convex minlps
    Options->GetIntegerValue("max_consecutive_failures",
        maxFailures,"bonmin.");
    Options->GetIntegerValue("max_consecutive_infeasible",
        maxInfeasible,"bonmin.");
    Options->GetEnumValue("nlp_failure_behavior",failureBehavior,"bonmin.");

    // Hybrid options
    Options->GetEnumValue("oa_cuts_scope", oaCutsGlobal,"bonmin.");
    Options->GetEnumValue("add_only_violated_oa", addOnlyViolatedOa,"bonmin.");
    Options->GetIntegerValue("nlp_solve_frequency",nlpSolveFrequency,"bonmin.");
    Options->GetIntegerValue("filmint_ecp_cuts",filmintCutsFrequency, "bonmin.");
    Options->GetIntegerValue("number_ecp_rounds", numEcpRounds,"bonmin.");
    Options->GetNumericValue("oa_dec_time_limit",oaDecMaxTime,"bonmin.");
    Options->GetIntegerValue("Gomory_cuts", migFreq,"bonmin.");
    Options->GetIntegerValue("probing_cuts",probFreq,"bonmin.");
    Options->GetIntegerValue("mir_cuts",mirFreq,"bonmin.");
    Options->GetIntegerValue("cover_cuts",coverFreq,"bonmin.");

    // milp subsolver options
    Options->GetEnumValue("milp_subsolver",milpSubSolver,"bonmin.");
    Options->GetEnumValue("nodeselect_stra",milpSubSolver_nodeSelection,"milp_sub.");
    Options->GetIntegerValue("number_strong_branch",milpSubSolver_numberStrong,"milp_sub.");
    Options->GetIntegerValue("number_before_trust", milpSubSolver_minReliability,"milp_sub.");
    Options->GetIntegerValue("Gomory_cuts", milpSubSolver_migFreq,"milp_sub.");
    Options->GetIntegerValue("probing_cuts",milpSubSolver_probFreq,"milp_sub.");
    Options->GetIntegerValue("mir_cuts",milpSubSolver_mirFreq,"milp_sub.");
    Options->GetIntegerValue("cover_cuts",milpSubSolver_coverFreq,"milp_sub.");

    //Preset default for algorithm
    if (algo==0)//B-BB
    {
    
    }
    else if (algo==1)//B-OA
    {
      oaDecMaxTime = DBL_MAX;
      nlpSolveFrequency = 0;
      bbLogLevel = 0;
    }
    else if (algo==2) {
      oaDecMaxTime = 0;
      nlpSolveFrequency = 0;
      migFreq = 0;
      probFreq = 0;
      mirFreq = 0;
      coverFreq = 0;
    }
    else if (algo==3)//Nothing to do
    {
    }
    
    // Set branching strategy
    if(varSelection == OsiTMINLPInterface::MOST_FRACTIONAL){
      minReliability = 0;
      numberStrong = 0;
    }
    else if(varSelection == OsiTMINLPInterface::STRONG_BRANCHING){
      minReliability = 0;
    }
    else if(varSelection == OsiTMINLPInterface::RELIABILITY_BRANCHING){
      minReliability = 10000;
    } 
    return success;
  }
}
