// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
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
    success &= Options->GetIntegerValue("bb_log_level",bbLogLevel,"bonmin.");
    success &= Options->GetIntegerValue("bb_log_interval",logInterval,"bonmin.");
    success &= Options->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    success &= Options->GetIntegerValue("milp_log_level",milpLogLevel,"bonmin.");
    success &= Options->GetIntegerValue("oa_log_level",oaLogLevel,"bonmin.");
    success &= Options->GetNumericValue("oa_log_frequency",oaLogFrequency,"bonmin.");
    success &= Options->GetIntegerValue("nlp_log_level",nlpLogLevel,"bonmin.");
    //General options
    success &= Options->GetEnumValue("algorithm",algo,"bonmin.");
    success &= Options->GetNumericValue("time_limit", maxTime, "bonmin.");
    success &= Options->GetIntegerValue("node_limit",maxNodes,"bonmin.");
    success &= Options->GetNumericValue("integer_tolerance",intTol,"bonmin.");
    success &= Options->GetNumericValue("allowable_gap",allowableGap,"bonmin.");
    success &= Options->GetNumericValue("allowable_fraction_gap",allowableFractionGap,"bonmin.");
    success &= Options->GetNumericValue("cutoff_decr",cutoffDecr,"bonmin.");
    success &= Options->GetNumericValue("cutoff",cutoff,"bonmin.");

    // Branch & bound setting
    success &= Options->GetEnumValue("nodeselect_stra",nodeSelection,"bonmin.");
    success &= Options->GetEnumValue("varselect_stra",varSelection,"bonmin.");
    success &= Options->GetIntegerValue("number_strong_branch",numberStrong,"bonmin.");
    success &= Options->GetIntegerValue("number_before_trust", minReliability,"bonmin.");
    success &= Options->GetIntegerValue("number_ecp_rounds", numEcpRounds,"bonmin.");

    success &=  Options->GetEnumValue("sos_constraints",disableSos,"bonmin.");
    // Robustness and non convex minlps
    success &= Options->GetIntegerValue("max_consecutive_failures",
        maxFailures,"bonmin.");
    success &= Options->GetIntegerValue("max_consecutive_infeasible",
        maxInfeasible,"bonmin.");
    success &= Options->GetEnumValue("nlp_failure_behavior",failureBehavior,".bonmin");

    // Hybrid options
    success &= Options->GetIntegerValue("nlp_solve_frequency",nlpSolveFrequency,"bonmin.");
success &= Options->GetIntegerValue("filmint_ecp_cuts",filmintCutsFrequency, "bonmin.");
    success &= Options->GetNumericValue("oa_dec_time_limit",oaDecMaxTime,"bonmin.");
    success &= Options->GetIntegerValue("Gomory_cuts", migFreq,"bonmin.");
    success &= Options->GetIntegerValue("probing_cuts",probFreq,"bonmin.");
    success &= Options->GetIntegerValue("mir_cuts",mirFreq,"bonmin.");
    success &= Options->GetIntegerValue("cover_cuts",coverFreq,"bonmin.");

    // milp subsolver options
    success &= Options->GetEnumValue("milp_subsolver",milpSubSolver,"bonmin.");
    success &= Options->GetEnumValue("nodeselect_stra",milpSubSolver_nodeSelection,"milp_sub.");
    success &= Options->GetIntegerValue("number_strong_branch",milpSubSolver_numberStrong,"milp_sub.");
    success &= Options->GetIntegerValue("number_before_trust", milpSubSolver_minReliability,"milp_sub.");
    success &= Options->GetIntegerValue("Gomory_cuts", milpSubSolver_migFreq,"milp_sub.");
    success &= Options->GetIntegerValue("probing_cuts",milpSubSolver_probFreq,"milp_sub.");
    success &= Options->GetIntegerValue("mir_cuts",milpSubSolver_mirFreq,"milp_sub.");
    success &= Options->GetIntegerValue("cover_cuts",milpSubSolver_coverFreq,"milp_sub.");

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
