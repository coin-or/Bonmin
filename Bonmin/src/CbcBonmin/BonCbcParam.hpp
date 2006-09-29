// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#ifndef BonminCbcParam_H
#define BonminCbcParam_H

#include "BonOsiTMINLPInterface.hpp"
namespace Bonmin
{
  class BonminCbcParam
  {
  public:
    /** Algorithm type
        <ul>
        <li> 0 "B-BB"
        <li> 1 "B-OA"
        <li> 2 "B-QG"
        <li> 3 "B-Hyb"
        </ul>
    */
    int algo;
    /** log level for the branch-and-bound */
    int bbLogLevel;
    /** Display information every logIntervval nodes.*/
    int logInterval;
    /** log level for the continuous subsolver */
    int lpLogLevel;
    /** log level for milp sub-solver in OA. */
    int milpLogLevel;
    /** log level for OA decomposition */
    int oaLogLevel;
    /** log frequency for OA */
    double oaLogFrequency;
    /** log level for the nlp subsolver interface (different from ipopt
        log and log level which should be set with print_level).*/
    int nlpLogLevel;
    /** Max number of failures in a branch.*/
    int maxFailures;
    /** Behavior of the algorithm in the case of a failure.*/
    int failureBehavior;
    /** Max number of consecutive infeasible problem in a branch
        before fathoming.*/
    int maxInfeasible;
    /** Amount by which cutoff is incremented */
    double cutoffDecr;
    /** cutoff value */
    double cutoff;
    /** Stop if absolute gap is less than :*/
    double allowableGap;
    /** Stop if relative gap is less than :*/
    double allowableFractionGap;
    /** Node selection strategy :
        <ul>
        <li> 0: best boud,
        <li> 1: DFS,
        <li> 2: BFS,
        <li> 3: dynamic (see
        <a href="http://www.coin-or.org/Doxygen/Cbc/class_cbc_branch_dynamic_decision.html">
        CbcBranchActual.hpp </a>
        </ul>
    */
    int nodeSelection;
    /** Number of candidates for strong branching.*/
    int numberStrong;
    /** Minimum reliability before trust pseudo-costs.*/
    int minReliability;
    /** Global time limit. */
    double maxTime;
    /** Global node limit.*/
    int maxNodes;
    /** Integer tolerance.*/
    double intTol;
    /** Conssider or not SOS constraints.*/
    int disableSos;
    /** frequency to solve nlp's in B-Hyb.*/
    int nlpSolveFrequency;
    /** Max OA decomposition time in B-Hyb.*/
    double oaDecMaxTime;
    /** milp subsolver:
        <ul>
        <li> 0 Cbc with defaults,
        <li> 1 Cbc with passed parameters,
        <li> 2 Cplex.
        </ul>
    */
    int milpSubSolver;
    /** Mig cuts generation frequency.*/
    int migFreq;
    /** Probing cuts generation frequency.*/
    int probFreq;
    /** Mir cuts generation frequency.*/
    int mirFreq;
    /** Cover cuts generation frequency.*/
    int coverFreq;

    /** (only set if milpSubSolver is 1) milpsubsolver
        Node selection strategy :
        <ul>
        <li> 0: best boud,
        <li> 1: DFS,
        <li> 2: BFS,
        <li> 3: dynamic (see
        <a href="http://www.coin-or.org/Doxygen/Cbc/class_cbc_branch_dynamic_decision.html">
        CbcBranchActual.hpp </a>
        </ul>
    */
    int milpSubSolver_nodeSelection;
    /** (only set if milpSubSolver is 1) milpsubsolver
        Number of candidates for strong branching.*/
    int milpSubSolver_numberStrong;
    /** (only set if milpSubSolver is 1) milpsubsolver
        Minimum reliability before trust pseudo-costs.*/
    int milpSubSolver_minReliability;
    /** (only set if milpSubSolver is 1) milpsubsolver Mig cuts generation frequency.*/
    int milpSubSolver_migFreq;
    /** (only set if milpSubSolver is 1) milpsubsolver Probing cuts generation frequency.*/
    int milpSubSolver_probFreq;
    /** (only set if milpSubSolver is 1) milpsubsolver Mir cuts generation frequency.*/
    int milpSubSolver_mirFreq;
    /** (only set if milpSubSolver is 1) milpsubsolver Cover cuts generation frequency.*/
    int milpSubSolver_coverFreq;

  public:
    /** Empty constructor. */
    BonminCbcParam()
    {}
    /** Destructor.*/
    ~BonminCbcParam()
    {}
    ///Process parameter file and extract MIP options.
    bool extractParams(OsiTMINLPInterface * solver);
    ///operator() will extractParameters from IpoptInterface.
    bool operator()(OsiTMINLPInterface * solver)
    {
      return extractParams(solver);
    }
  };
} //end namespace Bonmin
#endif

