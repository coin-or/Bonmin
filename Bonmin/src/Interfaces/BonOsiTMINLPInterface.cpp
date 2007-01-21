// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004

#include "BonminConfig.h"

#include "BonOsiTMINLPInterface.hpp"
#include "BonColReader.hpp"
#include "CoinTime.hpp"
#include <string>
#include <sstream>

#include "Ipopt/BonIpoptSolver.hpp"
#ifdef COIN_HAS_FILTERSQP
#include "Filter/BonFilterSolver.hpp"
#endif

using namespace Ipopt;


namespace Bonmin {
///Register options
static void
register_general_options
(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("bonmin output options");

  roptions->AddBoundedIntegerOption("bb_log_level",
      "specify main branch-and-bound log level.",
      0,5,1,
      "Set the level of output of the branch-and-bound : "
      "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                   );

  roptions->AddLowerBoundedIntegerOption("bb_log_interval",
      "Interval at which node level output is printed.",
      0,100,
      "Set the interval (in terms of number of nodes) at which "
      "a log on node resolutions (consisting of lower and upper bounds) is given.");

  roptions->AddBoundedIntegerOption("lp_log_level",
      "specify LP log level.",
      0,4,0,
      "Set the level of output of the linear programming sub-solver in B-Hyb or B-QG : "
      "0 - none, 1 - minimal, 2 - normal low, 3 - normal high, 4 - verbose"
                                   );

  roptions->AddBoundedIntegerOption("milp_log_level",
      "specify MILP subsolver log level.",
      0,3,0,
      "Set the level of output of the MILP subsolver in OA : "
      "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                   );

  roptions->AddBoundedIntegerOption("oa_log_level",
      "specify OA iterations log level.",
      0,2,1,
      "Set the level of output of OA decomposition solver : "
      "0 - none, 1 - normal, 2 - verbose"
                                   );

  roptions->AddLowerBoundedNumberOption("oa_log_frequency",
      "display an update on lower and upper bounds in OA every n seconds",
      0.,1.,100.,
      "");

  roptions->AddBoundedIntegerOption("nlp_log_level",
      "specify NLP solver interface log level (independent from ipopt print_level).",
      0,2,1,
      "Set the level of output of the OsiTMINLPInterface : "
      "0 - none, 1 - normal, 2 - verbose"
                                   );

  roptions->SetRegisteringCategory("bonmin branch-and-bound options");

  roptions->AddStringOption2("nlp_solver",
      "Choice of the solver for continuous nlp's",
      "Ipopt",
      "Ipopt", "Interior Point OPTimizer (https://projects.coin-or.org/Ipopt)",
      "filterSQP", "Sequential quadratic programming trust region algorithm (http://www-unix.mcs.anl.gov/~leyffer/solvers.html)",
       "");

  roptions->AddStringOption4("algorithm",
      "Choice of the algorithm.",
      "B-Hyb",
      "B-BB","simple branch-and-bound algorithm,",
      "B-OA","OA Decomposition algorithm,",
      "B-QG","Quesada and Grossmann branch-and-cut algorithm,",
      "B-Hyb","hybrid outer approximation based branch-and-cut.",
      "This will preset default values for most options of bonmin but depending on which algorithm "
      "some of these can be changed.");
  roptions->AddLowerBoundedNumberOption("time_limit",
      "Set the global maximum computation time (in secs) for the algorithm.",
      0.,1,1e10,
      "");

  roptions->AddLowerBoundedIntegerOption("node_limit",
      "Set the maximum number of nodes explored in the branch-and-bound search.",
      0,INT_MAX,
      "");

  roptions->AddBoundedNumberOption("integer_tolerance",
      "Set integer tolerance.",
      0.,1,.5,1,1e-06,
      "Any number within that value of an integer is considered integer.");

  roptions->AddBoundedNumberOption("allowable_gap",
      "Specify the value of absolute gap under which the algorithm stops.",
      -1.e20,0,1.e20,0,0.,
      "Stop the tree search when the gap between the objective value of the best known solution"
      " and the best bound on the objective of any solution is less than this.");

  roptions->AddBoundedNumberOption("allowable_fraction_gap",
      "Specify the value of relative gap under which the algorithm stops.",
      -1.e20,0,1.e20,0,0.0,
      "Stop the tree search when the gap between the objective value of the best known solution "
      "and the best bound on the objective of any solution is less than this "
      "fraction of the absolute value of the best known solution value.");

  roptions->AddBoundedNumberOption("cutoff",
      "Specify cutoff value.",
      -1e100,0,1e100,0,1e100,
      "cutoff should be the value of a feasible solution known by the user "
      "(if any). The algorithm will only look for solutions better than cutoof.");


  roptions->AddBoundedNumberOption("cutoff_decr",
      "Specify cutoff decrement.",
      -1e10,0,1e10,0,1e-05,
      "Specify the amount by which cutoff is decremented below "
      "a new best upper-bound"
      " (usually a small postive value but in non-convex problems it may be a negative value).");

  roptions->AddStringOption4("nodeselect_stra",
      "Choose the node selection strategy.",
      "best-bound",
      "best-bound", "choose node whith the smallest bound,",
      "depth-first", "Perform depth first search,",
      "breadth-first", "Perform breadth first search,",
      "dynamic", "Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 "
      "integer feasible solutions have been found).",
      "Choose the strategy for selecting the next node to be processed.");

  roptions->AddStringOption9("varselect_stra",
      "Chooses variable selection strategy",
      "strong-branching",
      "most-fractional", "Choose most fractional variable",
      "strong-branching", "Perform strong branching",
      "reliability-branching", "Use reliability branching",
      "curvature-estimator", "Use curvature estimation to select branching variable",
      "qp-strong-branching", "Perform strong branching with QP approximation",
      "lp-strong-branching", "Perform strong branching with LP approximation",
      "nlp-strong-branching", "Perform strong branching with NLP approximation",
      "osi-simple", "Osi method to do simple branching",
      "osi-strong", "Osi method to do strong branching","");

  roptions->AddLowerBoundedIntegerOption("number_strong_branch",
      "Choose the maximum number of variables considered for strong branching.",
      0,20,
      "Set the number of variables on which to do strong branching.");

  roptions->AddLowerBoundedIntegerOption
  ("number_before_trust",
   "Set the number of branches on a variable before its pseudo costs are to be believed "
   "in dynamic strong branching.",
   0,8,
   "A value of 0 disables dynamic strong branching.");

  roptions->AddLowerBoundedIntegerOption
  ("number_ecp_rounds",
   "Set the number of rounds of ecp in strong branching.",
   0,0,
   "");

  roptions->AddStringOption3("warm_start",
      "Select the warm start method",
      "optimum",
      "none","No warm start",
      "optimum","Warm start with direct parent optimum",
      "interior_point","Warm start with an interior point of direct parent",
      "This will affect the function getWarmStart(), and as a consequence the wam starting in the various algorithms.");

  roptions->AddStringOption2("sos_constraints",
			     "Wether or not to activate SOS constraints.",
			     "enable",
			     "enable","",
			     "disable","",
			     "(only type 1 SOS are supported at the moment)");

  roptions->SetRegisteringCategory("bonmin options for robustness");

  roptions->AddLowerBoundedNumberOption("max_random_point_radius",
      "Set max value r for coordinate of a random point.",
      0.,1,1e5,
      "When picking a random point coordinate i will be in the interval [min(max(l,-r),u-r), max(min(u,r),l+r)] "
      "(where l is the lower bound for the variable and u is its upper bound)");

  roptions->AddLowerBoundedIntegerOption
  ("max_consecutive_failures",
   "(temporarily removed) Number $n$ of consecutive unsolved problems before aborting a branch of the tree.",
   0,10,
   "When $n > 0$, continue exploring a branch of the tree until $n$ "
   "consecutive problems in the branch are unsolved (we call unsolved a problem for which Ipopt can not "
   "guarantee optimality within the specified tolerances).");

  roptions->AddLowerBoundedIntegerOption
  ("num_iterations_suspect",
   "Number of iterations over which a node is considered \"suspect\" (for debugging purposes only, see detailed documentation).",
   -1,-1,
   "When the number of iterations to solve a node is above this number, the subproblem at this"
   " node is considered to be suspect and it will be outputed in a file (set to -1 to deactivate this).");

  roptions->AddStringOption2("nlp_failure_behavior",
      "Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not "
      "able to guarantee optimality within the specified tolerances).",
      "stop",
      "stop", "Stop when failure happens.",
      "fathom", "Continue when failure happens.",
      "If set to \"fathom\", the algorithm will fathom the node when Ipopt fails to find a solution to the nlp "
      "at that node whithin the specified tolerances. "
      "The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.");

  roptions->AddLowerBoundedIntegerOption("num_retry_unsolved_random_point",
      "Number $k$ of times that the algorithm will try to resolve an unsolved NLP with a random starting point "
      "(we call unsolved an NLP for which Ipopt is not "
      "able to guarantee optimality within the specified tolerances).",
      0,0,
      "When Ipopt fails to solve a continuous NLP sub-problem, if $k > 0$, the algorithm will "
      "try again to solve the failed NLP with $k$ new randomly chosen starting points "
      " or until the problem is solved with success.");


  roptions->SetRegisteringCategory("bonmin options for non-convex problems");
  roptions->AddLowerBoundedIntegerOption("max_consecutive_infeasible",
      "Number of consecutive infeasible subproblems before aborting a"
      " branch.",
      0,0,
      "Will continue exploring a branch of the tree until \"max_consecutive_infeasible\""
      "consecutive problems are infeasibles by the NLP sub-solver.");

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_root",
      "Number $k$ of tries to resolve the root node with different starting points.",
      0,0,
      "The algorithm will solve the root node with $k$ random starting points"
      " and will keep the best local optimum found.");

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_node",
      "Number $k$ of tries to resolve a node (other than the root) of the tree with different starting point.",
      0,0,
      "The algorithm will solve all the nodes with $k$ different random starting points "
      "and will keep the best local optimum found.");


}

static void register_OA_options
(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("bonmin options : B-Hyb specific options");

  roptions->AddLowerBoundedIntegerOption("nlp_solve_frequency",
      "Specify the frequency (in terms of nodes) at which NLP relaxations are solved in B-Hyb.",
      0,10,
      "A frequency of 0 amounts to to never solve the NLP relaxation.");

  roptions->AddLowerBoundedIntegerOption("filmint_ecp_cuts",
      "Specify the frequency (in terms of nodes) at which some a la filmint ecp cuts are generated.",
      0,0,
      "A frequency of 0 amounts to to never solve the NLP relaxation.");

  roptions->AddLowerBoundedNumberOption("oa_dec_time_limit",
      "Specify the maximum number of seconds spent overall in OA decomposition iterations.",
      0.,0,30.,
      "");
  roptions->AddLowerBoundedNumberOption("tiny_element","Value for tiny element in OA cut",
      -0.,0,1e-08,
      "We will remove \"cleanly\" (by relaxing cut) an element lower"
      " than this.");

  roptions->AddLowerBoundedNumberOption("very_tiny_element","Value for very tiny element in OA cut",
      -0.,0,1e-17,
      "Algorithm will take the risk of neglecting an element lower"
      " than this.");

  roptions->AddLowerBoundedIntegerOption("Gomory_cuts",
      "Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.",
      -100,-5,
      "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
      "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
      "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
  roptions->AddLowerBoundedIntegerOption("probing_cuts",
      "Frequency (in terms of nodes) for generating probing cuts in branch-and-cut",
      -100,-5,
      "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
      "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
      "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");

  roptions->AddLowerBoundedIntegerOption("cover_cuts",
      "Frequency (in terms of nodes) for generating cover cuts in branch-and-cut",
      -100,-5,
      "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
      "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
      "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");

  roptions->AddLowerBoundedIntegerOption("mir_cuts",
      "Frequency (in terms of nodes) for generating MIR cuts in branch-and-cut",
      -100,-5,
      "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
      "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
      "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");

}

static void register_milp_sub_solver_options
(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("bonmin options : Options for milp subsolver in OA decomposition");
  roptions->AddStringOption3("milp_subsolver",
      "Choose the subsolver to solve MILPs sub-problems in OA decompositions.",
      "Cbc_D",
      "Cbc_D","Coin Branch and Cut with its default",
      "Cbc_Par", "Coin Branch and Cut with passed parameters",
      "Cplex","Ilog Cplex",
      " To use Cplex, a valid license is required and you should have compiled OsiCpx in COIN-OR  (see Osi documentation).");
}

///Register options
void
OsiTMINLPInterface::register_ALL_options
(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions)
{
  // We try to register the options - if those have been registered
  // already, we catch the exception and don't need to do it again
  try {
    register_general_options(roptions);
    register_OA_options(roptions);
    register_milp_sub_solver_options(roptions);
    //Register options for all possible solvers (besides the one used
    IpoptSolver * ipopt = dynamic_cast<IpoptSolver *> (GetRawPtr(app_));
#ifdef COIN_HAS_FILTERSQP
    FilterSolver * filter = dynamic_cast<FilterSolver *> (GetRawPtr(app_));
#endif
    if(!ipopt) {
      ipopt = new IpoptSolver;
      ipopt->RegisterOptions(app_->RegOptions());
      delete ipopt;
      ipopt = NULL;
    }

#ifdef COIN_HAS_FILTERSQP
    if(!filter) {
      filter = new FilterSolver;
      filter->RegisterOptions(app_->RegOptions());
      delete filter;
      filter = NULL;
    }
#endif		
  }   
  catch(RegisteredOptions::OPTION_ALREADY_REGISTERED) {
    // skipping
  }
}


OsiTMINLPInterface::Messages::Messages
():CoinMessages((int)OSITMINLPINTERFACE_DUMMY_END)
{
  strcpy(source_ ,"NLP");
  addMessage(SOLUTION_FOUND, CoinOneMessage
      (1,2,"After %d tries found a solution of %g "
       "(previous best %g)."));
  addMessage(INFEASIBLE_SOLUTION_FOUND, CoinOneMessage
      (2,2,"After %d tries found an solution of %g "
       "infeasible problem."));

  addMessage(UNSOLVED_PROBLEM_FOUND, CoinOneMessage
      (3,2,"After %d tries found an solution of %g "
       "unsolved problem."));
  addMessage(WARN_SUCCESS_WS,CoinOneMessage
      (3009,2,
       "Problem not solved with warm start but "
       "solved without"));

  addMessage(WARNING_RESOLVING,CoinOneMessage
      (3010,2,
       "Trying to resolve NLP with different starting "
       "point (%d attempts)."));
  addMessage(WARN_SUCCESS_RANDOM, CoinOneMessage
      (3011,1,
       "Problem initially not solved but solved with "
       "a random starting point (success on %d attempt)"));
  addMessage(WARN_CONTINUING_ON_FAILURE, CoinOneMessage
      (3017,1,
       "Warning : continuing branching, while there are "
       "unrecovered failures at nodes"));

  addMessage(SUSPECT_PROBLEM, CoinOneMessage
      (4,2,"NLP number %d is suspect (see bounds and start file)"));
  addMessage(IPOPT_SUMMARY, CoinOneMessage
      (6,2,"Ipopt return (for %s): status %2d, iter count %4d, time %g"));
  addMessage(BETTER_SOL, CoinOneMessage
      (7,2,"Solution of value %g found on %d'th attempt"));

  addMessage(LOG_HEAD, CoinOneMessage
      (8,1,"\n          "
       "    Num      Status      Obj             It       time"));
  addMessage(LOG_FIRST_LINE, CoinOneMessage
      (9,1,
       "    %-8d %-11s %-14g %-8d %-3g"));
  addMessage(LOG_LINE, CoinOneMessage
      (10,1,
       " %c  r%-7d %-11s %-14g %-8d %-3g"));

  addMessage(WARN_RESOLVE_BEFORE_INITIAL_SOLVE, CoinOneMessage
      (3012,1,"resolve called before any call to initialSolve"
       " can not use warm starts."));

}
bool OsiTMINLPInterface::hasPrintedOptions=0;

////////////////////////////////////////////////////////////////////
// Constructors and desctructors                                  //
////////////////////////////////////////////////////////////////////
/// Default Constructor
OsiTMINLPInterface::OsiTMINLPInterface():
    OsiSolverInterface(),
    tminlp_(NULL),
    problem_(NULL),
    app_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    reducedCosts_(NULL),
    OsiDualObjectiveLimit_(1e200),
//   optimization_status_(HasNotBeenOptimized),
    varNames_(NULL),
    hasVarNamesFile_(true),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(1e08),
    pushValue_(1e-02),
    numRetryInitial_(-1),
    numRetryResolve_(-1),
    numRetryUnsolved_(1),
    messages_(),
    pretendFailIsInfeasible_(0),
    hasContinuedAfterNlpFailure_(false),
    numIterationSuspect_(-1),
    hasBeenOptimized_(false),
    obj_(NULL),
    feasibilityProblem_(NULL),
    jRow_(NULL),
    jCol_(NULL),
    jValues_(NULL),
    nnz_jac(0),
    constTypes_(NULL),
//    constTypesNum_(NULL),
    nLinear_(0),
    nNonLinear_(0),
    tiny_(1e-8),
    veryTiny_(1e-20),
    infty_(1e100),
    expose_warm_start_(false)
{}

/** Constructor with given TNLPSolver and TMINLP */
OsiTMINLPInterface::OsiTMINLPInterface (Ipopt::SmartPtr<Bonmin::TNLPSolver> app):
    OsiSolverInterface(),
    tminlp_(NULL),
    problem_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    reducedCosts_(NULL),
    OsiDualObjectiveLimit_(1e200),
//    optimization_status_(HasNotBeenOptimized),
    varNames_(NULL),
    hasVarNamesFile_(true),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(1e08),
    pushValue_(1e-02),
    numRetryInitial_(-1),
    numRetryResolve_(-1),
    numRetryUnsolved_(1),
    messages_(),
    pretendFailIsInfeasible_(false),
    hasContinuedAfterNlpFailure_(false),
    numIterationSuspect_(-1),
    hasBeenOptimized_(false),
    obj_(NULL),
    feasibilityProblem_(NULL),
    jRow_(NULL),
    jCol_(NULL),
    jValues_(NULL),
    nnz_jac(0),
    constTypes_(NULL),
//    constTypesNum_(NULL),
    nLinear_(0),
    nNonLinear_(0),
    tiny_(1e-08),
    veryTiny_(1e-17),
    infty_(1e100),
    expose_warm_start_(false)
{
  app_ = app->clone();

  Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = app_->RegOptions();
  register_ALL_options(roptions);
  app_->Initialize("bonmin.opt");
  extractInterfaceParams();

}

/** Constructor with given IpSolver and TMINLP */
OsiTMINLPInterface::OsiTMINLPInterface (Ipopt::SmartPtr<Bonmin::TMINLP> tminlp,
                                        Ipopt::SmartPtr<Bonmin::TNLPSolver> app):
    OsiSolverInterface(),
    tminlp_(tminlp),
    problem_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    reducedCosts_(NULL),
    OsiDualObjectiveLimit_(1e200),
//    optimization_status_(HasNotBeenOptimized),
    varNames_(NULL),
    hasVarNamesFile_(true),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(1e08),
    pushValue_(1e-02),
    numRetryInitial_(-1),
    numRetryResolve_(-1),
    numRetryUnsolved_(1),
    messages_(),
    pretendFailIsInfeasible_(false),
    hasContinuedAfterNlpFailure_(false),
    numIterationSuspect_(-1),
    hasBeenOptimized_(false),
    obj_(NULL),
    feasibilityProblem_(NULL),
    jRow_(NULL),
    jCol_(NULL),
    jValues_(NULL),
    nnz_jac(0),
    constTypes_(NULL),
//    constTypesNum_(NULL),
    nLinear_(0),
    nNonLinear_(0),
    tiny_(1e-08),
    veryTiny_(1e-17),
    infty_(1e100),
    expose_warm_start_(false)
{
  allocateTMINLP(tminlp,app);
}
/** Facilitator to allocate a tminlp and an application. */
void
OsiTMINLPInterface::allocateTMINLP(Ipopt::SmartPtr<Bonmin::TMINLP> tminlp,
				   Ipopt::SmartPtr<Bonmin::TNLPSolver> app)
{
  assert(IsValid(tminlp));

  tminlp_ = tminlp;
  problem_ = new TMINLP2TNLP(tminlp_);
  app_ = app->clone();

  Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = app_->RegOptions();
  register_ALL_options(roptions);
  app_->Initialize("");
  extractInterfaceParams();
}




void
OsiTMINLPInterface::readOptionFile(const char * fileName)
{
  extractInterfaceParams();
}

/// Copy constructor
OsiTMINLPInterface::OsiTMINLPInterface (const OsiTMINLPInterface &source):
    OsiSolverInterface(source),
    tminlp_(source.tminlp_),
    problem_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    reducedCosts_(NULL),
    OsiDualObjectiveLimit_(source.OsiDualObjectiveLimit_),
    varNames_(NULL),
    hasVarNamesFile_(source.hasVarNamesFile_),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(source.maxRandomRadius_),
    pushValue_(source.pushValue_),
    numRetryInitial_(source.numRetryInitial_),
    numRetryResolve_(source.numRetryResolve_),
    numRetryUnsolved_(source.numRetryUnsolved_),
    messages_(),
    pretendFailIsInfeasible_(source.pretendFailIsInfeasible_),
    hasContinuedAfterNlpFailure_(source.hasContinuedAfterNlpFailure_),
    numIterationSuspect_(source.numIterationSuspect_),
    hasBeenOptimized_(source.hasBeenOptimized_),
    obj_(NULL),
    feasibilityProblem_(NULL),
    jRow_(NULL),
    jCol_(NULL),
    jValues_(NULL),
    nnz_jac(source.nnz_jac),
    constTypes_(NULL),
//    constTypesNum_(NULL),
    nLinear_(0),
    nNonLinear_(0),
    tiny_(source.tiny_),
    veryTiny_(source.veryTiny_),
    infty_(source.infty_),
    expose_warm_start_(source.expose_warm_start_)
{
  // Copy options from old application
  if(IsValid(source.tminlp_)) {
    problem_ = new Bonmin::TMINLP2TNLP(tminlp_);
    problem_->copyUserModification(*source.problem_);
    pretendFailIsInfeasible_ = source.pretendFailIsInfeasible_;

    // Copy options from old application
    app_ = source.app_->clone();
  }
  else {
    throw SimpleError("Don't know how to copy an empty IpoptInterface.",
        "copy constructor");
  }

  if(source.obj_) {
    obj_ = new double[source.getNumCols()];
    CoinCopyN(source.obj_, source.getNumCols(), obj_);
  }
  if(IsValid(source.tminlp_))
    feasibilityProblem_ = new Bonmin::TNLP2FPNLP
        (Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(problem_)));
  else
    throw SimpleError("Don't know how to copy an empty IpoptOAInterface.",
        "copy constructor");


  if(source.jValues_!=NULL && source.jRow_ != NULL && source.jCol_ != NULL && nnz_jac>0) {
    jValues_ = new double [nnz_jac];
    jCol_    = new Ipopt::Index [nnz_jac];
    jRow_    = new Ipopt::Index [nnz_jac];
    CoinCopyN(source.jValues_ , nnz_jac,jValues_ );
    CoinCopyN(source.jCol_    , nnz_jac,jCol_    );
    CoinCopyN(source.jRow_    , nnz_jac,jRow_    );

    if(source.constTypes_ != NULL) {
      constTypes_ = new Ipopt::TNLP::LinearityType [getNumRows()];
      CoinCopyN(source.constTypes_, getNumRows(), constTypes_);
    }
  }
  else if(nnz_jac > 0) {
    throw SimpleError("Arrays for storing jacobian are inconsistant.",
        "copy constructor");
  }
}

OsiSolverInterface * 
OsiTMINLPInterface::clone(bool copyData ) const
{
  if(copyData)
    return new OsiTMINLPInterface(*this);
  else return new OsiTMINLPInterface;
}

/// Assignment operator
OsiTMINLPInterface & OsiTMINLPInterface::operator=(const OsiTMINLPInterface& rhs)
{
  if(this!= &rhs) {
    OsiSolverInterface::operator=(rhs);
    OsiDualObjectiveLimit_ = rhs.OsiDualObjectiveLimit_;
    nCallOptimizeTNLP_ = rhs.nCallOptimizeTNLP_;
    totalNlpSolveTime_ = rhs.nCallOptimizeTNLP_;
    totalIterations_ = rhs.totalIterations_;
    maxRandomRadius_ = rhs.maxRandomRadius_;
    hasVarNamesFile_ = rhs.hasVarNamesFile_;
    pushValue_ = rhs.pushValue_;

    if(IsValid(rhs.tminlp_)) {

      tminlp_ = rhs.tminlp_;
      problem_ = new Bonmin::TMINLP2TNLP(tminlp_);
      app_ = rhs.app_->clone();

      feasibilityProblem_ = new Bonmin::TNLP2FPNLP
          (Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(problem_)));
      nnz_jac = rhs.nnz_jac;

      if(constTypes_ != NULL) {
        delete [] constTypes_;
        constTypes_ = NULL;
      }
      if(rhs.constTypes_ != NULL) {
        constTypes_ = new Ipopt::TNLP::LinearityType[getNumRows()];
        CoinCopyN(rhs.constTypes_, getNumRows(), constTypes_);
      }
/*
      if(constTypesNum_ != NULL) {
        delete [] constTypesNum_;
        constTypesNum_ = NULL;
      }
      if(rhs.constTypesNum_ != NULL) {
        constTypesNum_ = new int[getNumRows()];
        CoinCopyN(rhs.constTypesNum_, getNumRows(), constTypesNum_);
      }
*/
      if(rhs.jValues_!=NULL && rhs.jRow_ != NULL && rhs.jCol_ != NULL && nnz_jac>0) {
        jValues_ = new double [nnz_jac];
        jCol_    = new Ipopt::Index [nnz_jac];
        jRow_    = new Ipopt::Index [nnz_jac];
        CoinCopyN(rhs.jValues_ , nnz_jac,jValues_ );
        CoinCopyN(rhs.jCol_    , nnz_jac,jCol_    );
        CoinCopyN(rhs.jRow_    , nnz_jac,jRow_    );
      }
      else if(nnz_jac > 0) {
        throw CoinError("Arrays for storing jacobian are inconsistant.",
            "copy constructor",
            "IpoptOAInterface");
      }
      tiny_ = rhs.tiny_;
      veryTiny_ = rhs.veryTiny_;
      infty_ = rhs.infty_;
      expose_warm_start_ = rhs.expose_warm_start_;

    }
    else {
      tminlp_ =NULL;
      problem_ = NULL;
app_ = NULL;
      feasibilityProblem_ = NULL;
    }


    if(obj_) {
      delete [] obj_;
      obj_ = NULL;
    }
    if(rhs.obj_) {
      obj_ = new double[rhs.getNumCols()];
      CoinCopyN(rhs.obj_, rhs.getNumCols(), obj_);
    }

    delete [] varNames_;
    varNames_ = NULL;

    if(rhs.varNames_) {
      rhs.varNames_ = new std::string[getNumCols()];
      CoinCopyN(rhs.varNames_, getNumCols(), varNames_);
    }

    hasVarNamesFile_ = rhs.hasVarNamesFile_;

    nCallOptimizeTNLP_ = rhs.nCallOptimizeTNLP_;
    totalNlpSolveTime_ = rhs.totalNlpSolveTime_;
    totalIterations_ = rhs.totalIterations_;
    maxRandomRadius_ = rhs.maxRandomRadius_;
    pushValue_ = rhs.pushValue_;
    numRetryInitial_ = rhs.numRetryInitial_;
    numRetryResolve_ = rhs.numRetryResolve_;
    numRetryUnsolved_ = rhs.numRetryUnsolved_;
    pretendFailIsInfeasible_ = rhs.pretendFailIsInfeasible_;
    numIterationSuspect_ = rhs.numIterationSuspect_;

    hasBeenOptimized_ = rhs.hasBeenOptimized_;

    freeCachedData();
  }
  return *this;
}

Ipopt::SmartPtr<Ipopt::OptionsList> OsiTMINLPInterface::retrieve_options()
{
  if(!IsValid(app_)) {
    std::cout<<"Can not parse options when no IpApplication has been created"<<std::endl;
    return NULL;
  }
  else
    return app_->Options();
}

/// Destructor
OsiTMINLPInterface::~OsiTMINLPInterface ()
{
  freeCachedData();
  delete [] jRow_;
  delete [] jCol_;
  delete [] jValues_;
  delete [] constTypes_;
  delete [] obj_;
}

void
OsiTMINLPInterface::freeCachedColRim()
{
  if(varNames_!=NULL) {
    delete [] varNames_;
    varNames_ = NULL;
  }
  if(reducedCosts_!=NULL) {
    delete []  reducedCosts_;
    reducedCosts_ = NULL;
  }

}

void
OsiTMINLPInterface::freeCachedRowRim()
{
  if(rowsense_!=NULL) {
    delete []  rowsense_;
    rowsense_ = NULL;
  }
  if(rhs_!=NULL) {
    delete []  rhs_;
    rhs_ = NULL;
  }
  if(rowrange_!=NULL) {
    delete []  rowrange_;
    rowrange_ = NULL;
  }
  //   if(dualsol_!=NULL)
  //     {
  //       delete []  dualsol_; dualsol_ = NULL;
  //     }
}

void
OsiTMINLPInterface::freeCachedData()
{
  freeCachedColRim();
  freeCachedRowRim();
}
static const char * OPT_SYMB="OPT";
static const char * FAILED_SYMB="FAILED";
static const char * INFEAS_SYMB="INFEAS";

///////////////////////////////////////////////////////////////////
// WarmStart Information                                                                           //
///////////////////////////////////////////////////////////////////


void
OsiTMINLPInterface::resolveForCost(int numsolve)
{

  /** Save current optimum. */
  double * point = new double[getNumCols()*3+ getNumRows()];
  double bestBound = getObjValue();
  CoinCopyN(getColSolution(),
      getNumCols(), point);
  CoinCopyN(getRowPrice(),
      2*getNumCols()+ getNumRows(),
      &point[getNumCols()]);

  if(isProvenOptimal())
    messageHandler()->message(SOLUTION_FOUND,
        messages_)
    <<1<<getObjValue()<<bestBound
    <<CoinMessageEol;
  else
    messageHandler()->message(INFEASIBLE_SOLUTION_FOUND,
        messages_)
    <<1
    <<CoinMessageEol;
  for(int f = 0; f < numsolve ; f++) {
    messageHandler()->message(WARNING_RESOLVING,
        messages_)
    <<f+1<< CoinMessageEol ;
    randomStartingPoint();
    solveAndCheckErrors(0,0,"resolveForCost");


    const char * status=OPT_SYMB;
    ;
    char c=' ';
    if(isAbandoned()) {
      status=FAILED_SYMB;
    }
    else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;


    //Is solution better than previous
    if(isProvenOptimal() &&
        getObjValue()<bestBound) {
      c='*';
      messageHandler()->message(BETTER_SOL, messages_)<<getObjValue()<<f+1<< CoinMessageEol;
      CoinCopyN(getColSolution(),
          getNumCols(), point);
      CoinCopyN(getRowPrice(),
          2*getNumCols()+ getNumRows(),
          &point[getNumCols()]);
      bestBound = getObjValue();
    }

    messageHandler()->message(LOG_LINE, messages_)
    <<c<<f+1<<status<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    if(isProvenOptimal())
      messageHandler()->message(SOLUTION_FOUND,
          messages_)
      <<f+2<<getObjValue()<<bestBound
      <<CoinMessageEol;
    else if(!isAbandoned())
      messageHandler()->message(UNSOLVED_PROBLEM_FOUND,
          messages_)
      <<f+2
      <<CoinMessageEol;
    else
      messageHandler()->message(INFEASIBLE_SOLUTION_FOUND,
          messages_)
      <<f+2
      <<CoinMessageEol;
  }
  //  setColSolution(point);
  //  setRowPrice(&point[getNumCols()]);
  //  setWarmStartOptions();
  delete [] point;

  optimization_status_ = app_->ReOptimizeTNLP(GetRawPtr(problem_));
  hasBeenOptimized_ = true;
}

void
OsiTMINLPInterface::resolveForRobustness(int numsolve)
{
  //std::cerr<<"Resolving the problem for robustness"<<std::endl;
  //First remove warm start point and resolve
  app_->disableWarmStart();
  messageHandler()->message(WARNING_RESOLVING,
      messages_)
  <<1<< CoinMessageEol ;
  solveAndCheckErrors(0,0,"resolveForRobustness");


  const char * status=OPT_SYMB;
  ;
  char c='*';
  if(isAbandoned()) {
    status=FAILED_SYMB;
    c=' ';
  }
  else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;
  messageHandler()->message(LOG_LINE, messages_)
  <<c<<1<<status<<getObjValue()<<app_->IterationCount()<<
    app_->CPUTime()<<CoinMessageEol;


  if(!isAbandoned()) {
    messageHandler()->message(WARN_SUCCESS_WS,
        messages_)
    << CoinMessageEol ;
    return; //we won go on
  }

  //still unsolved try again with different random starting points
  for(int f = 0; f < numsolve ; f++) {
    messageHandler()->message(WARNING_RESOLVING,
        messages_)
    <<f+2<< CoinMessageEol ;

    randomStartingPoint();
    solveAndCheckErrors(0,0,"resolveForRobustness");


    messageHandler()->message(IPOPT_SUMMARY, messages_)
    <<"resolveForRobustness"<<optimization_status_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    const char * status=OPT_SYMB;
    ;
    char c='*';
    if(isAbandoned()) {
      status=FAILED_SYMB;
      c=' ';
    }
    else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;
    messageHandler()->message(LOG_LINE, messages_)
    <<c<<f+2<<status<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


    if(!isAbandoned()) {
      messageHandler()->message(WARN_SUCCESS_RANDOM,
          messages_)
      <<f+2
      << CoinMessageEol ;
      return; //we have found a solution and continue
    }
  }
  if(pretendFailIsInfeasible_) {
    if(pretendFailIsInfeasible_ == 1) {
      messageHandler()->message(WARN_CONTINUING_ON_FAILURE,
          messages_)
      <<CoinMessageEol;
      hasContinuedAfterNlpFailure_ = 1;
    }
    return;
  }
  else {
    std::string probName;
    getStrParam(OsiProbName,probName);
    throw newUnsolvedError(optimization_status_, problem_,
                           probName);
  }
}

////////////////////////////////////////////////////////////////////
// Problem information methods                                    //
////////////////////////////////////////////////////////////////////
/// Get number of columns
int OsiTMINLPInterface::getNumCols() const
{

  return problem_->num_variables();
}


/// Get number of rows
int
OsiTMINLPInterface::getNumRows() const
{
  return problem_->num_constraints();
}

const double *
OsiTMINLPInterface::getColLower() const
{
  return problem_->x_l();
}

const double *
OsiTMINLPInterface::getColUpper() const
{
  return problem_->x_u();
}

#if 0
void
OsiTMINLPInterface::readVarNames() const
{
  delete []varNames_;
  varNames_ = NULL;
  std::string probName;
  getStrParam(OsiProbName, probName);
  ColReader colRead(probName);
  if(colRead.readFile()) {
    varNames_ = new std::string[problem_->num_variables()];
    colRead.copyNames(varNames_,problem_->num_variables());
    hasVarNamesFile_ = true;
  }
  else
    hasVarNamesFile_ = false;
}

///get name of a variable
const std::string *
OsiTMINLPInterface::getVarNames() const
{
  if(varNames_ == NULL && hasVarNamesFile_ ) {
    readVarNames();
  }
  return varNames_;
}
#endif


void OsiTMINLPInterface::extractSenseRhsAndRange() const
{
  assert(rowsense_==NULL&&rhs_==NULL&&rowrange_==NULL);
  int numrows = problem_->num_constraints();
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  rowsense_ = new char [numrows];
  rhs_ = new double [numrows];
  rowrange_ = new double [numrows];
  for(int i = 0 ; i < numrows ; i++) {
    rowrange_[i]=0.;
    convertBoundToSense(rowLower[i], rowUpper[i], rowsense_[i], rhs_[i], rowrange_[i]);
  }
}

/** Get pointer to array[getNumRows()] of row constraint senses.
    <ul>
    <li>'L': <= constraint
    <li>'E': =  constraint
    <li>'G': >= constraint
    <li>'R': ranged constraint
    <li>'N': free constraint
    </ul>
*/
const char *
OsiTMINLPInterface::getRowSense() const
{
  if(rowsense_==NULL) {
    extractSenseRhsAndRange();
  }
  return rowsense_;
}

/** Get pointer to array[getNumRows()] of rows right-hand sides
    <ul>
    <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
    <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
    <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
    <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
    </ul>
*/
const double *
OsiTMINLPInterface::getRightHandSide() const
{
  if(rhs_==NULL) {
    extractSenseRhsAndRange();
  }
  return rhs_;
}

/** Get pointer to array[getNumRows()] of row ranges.
    <ul>
    <li> if rowsense()[i] == 'R' then
    rowrange()[i] == rowupper()[i] - rowlower()[i]
    <li> if rowsense()[i] != 'R' then
    rowrange()[i] is 0.0
    </ul>
*/
const double *
OsiTMINLPInterface::getRowRange() const
{
  if(rowrange_==NULL) {
    extractSenseRhsAndRange();
  }
  return rowrange_;
}

const double *
OsiTMINLPInterface::getRowLower() const
{
  return problem_->g_l();
}

const double *
OsiTMINLPInterface::getRowUpper() const
{
  return problem_->g_u();
}

/// Return true if column is continuous
bool
OsiTMINLPInterface::isContinuous(int colNumber) const
{
  return (problem_->var_types()[colNumber]==Bonmin::TMINLP::CONTINUOUS);
}

/// Return true if column is binary
bool
OsiTMINLPInterface::isBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==Bonmin::TMINLP::BINARY);
}

/** Return true if column is integer.
    Note: This function returns true if the the column
    is binary or a general integer.
*/
bool
OsiTMINLPInterface::isInteger(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==Bonmin::TMINLP::BINARY)||
      (problem_->var_types()[colNumber]==Bonmin::TMINLP::INTEGER));
}

/// Return true if column is general integer
bool
OsiTMINLPInterface::isIntegerNonBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==Bonmin::TMINLP::INTEGER);
}
/// Return true if column is binary and not fixed at either bound
bool
OsiTMINLPInterface::isFreeBinary(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==Bonmin::TMINLP::BINARY)
      &&((getColUpper()[colNumber]-getColLower()[colNumber]) > 1 - 1e-09));
}

/// Get solver's value for infinity
double
OsiTMINLPInterface::getInfinity() const
{
  return DBL_MAX;
}

/// Get pointer to array[getNumCols()] of primal solution vector
const double *
OsiTMINLPInterface::getColSolution() const
{
  if(hasBeenOptimized_)
    return problem_->x_sol();
  else
    return problem_->x_init();
}

/// Get pointer to array[getNumRows()] of dual prices
const double *
OsiTMINLPInterface::getRowPrice() const
{
  if(hasBeenOptimized_)
    return problem_->duals_sol();
  else
    return problem_->duals_init();
}

/// Get a pointer to array[getNumCols()] of reduced costs
const double *
OsiTMINLPInterface::getReducedCost() const
{
  std::cerr<<"WARNING : trying to access reduced cost in Ipopt always retrun 0"<<std::endl;
  if(reducedCosts_==NULL) {
    reducedCosts_ = new double [getNumCols()];
    CoinFillN(reducedCosts_,getNumCols(),0.);
  }
  return reducedCosts_;
}

/** Get pointer to array[getNumRows()] of row activity levels (constraint
    matrix times the solution vector */
const double *
OsiTMINLPInterface::getRowActivity() const
{
  return problem_->g_sol();
}

/** Get how many iterations it took to solve the problem (whatever
    "iteration" mean to the solver.
*/
int
OsiTMINLPInterface::getIterationCount() const
{
    return app_->IterationCount();
}


/** Set a single column lower bound.
    Use -getInfinity() for -infinity. */
void
OsiTMINLPInterface::setColLower( int elementIndex, double elementValue )
{
  //  if(fabs(problem_->x_l()[elementIndex]-elementValue)>1e-06)
  problem_->SetVariableLowerBound(elementIndex,elementValue);
  hasBeenOptimized_ = false;
}

/** Set a single column upper bound.
    Use getInfinity() for infinity. */
void
OsiTMINLPInterface::setColUpper( int elementIndex, double elementValue )
{
  //  if(fabs(problem_->x_u()[elementIndex]-elementValue)>1e-06)
  problem_->SetVariableUpperBound(elementIndex,elementValue);
  hasBeenOptimized_ = false;
}

/** Set the lower bounds for all columns
    Use -getInfinity() for -infinity. */
void
OsiTMINLPInterface::setColLower( const double* array )
{
  problem_->SetVariablesLowerBounds(problem_->num_variables(),
                                  array);
  hasBeenOptimized_ = false;
}

/** Set Set the upper bounds for all columns
    Use getInfinity() for infinity. */
void
OsiTMINLPInterface::setColUpper( const double* array )
{
  problem_->SetVariablesUpperBounds(problem_->num_variables(), 
                                  array);
  hasBeenOptimized_ = false;
}

/** Set a single row lower bound.
    Use -getInfinity() for -infinity. */
void
OsiTMINLPInterface::setRowLower( int elementIndex, double elementValue )
{
  throw SimpleError("Not implemented yet but should be if necessary.",
      "setRowLower");
  hasBeenOptimized_ = false;
}

/** Set a single row upper bound.
    Use getInfinity() for infinity. */
void
OsiTMINLPInterface::setRowUpper( int elementIndex, double elementValue )
{
  throw SimpleError("Not implemented yet but should be if necessary.",
      "setRowUpper");
  hasBeenOptimized_ = false;
}

/** Set the type of a single row */
void
OsiTMINLPInterface::setRowType(int index, char sense, double rightHandSide,
    double range)
{
  throw SimpleError("Not implemented yet but should be if necessary.",
      "setRowType");
  hasBeenOptimized_ = false;
}


/// Set the objective function sense.
/// (1 for min (default), -1 for max)
void
OsiTMINLPInterface::setObjSense(double s)
{
  throw SimpleError("Can not change objective sense of an Ipopt problem.",
      "setObjSense");
  hasBeenOptimized_ = false;
}

/** Set the primal solution variable values

colsol[getNumCols()] is an array of values for the primal variables.
These values are copied to memory owned by the solver interface object
or the solver.  They will be returned as the result of getColSolution()
until changed by another call to setColSolution() or by a call to any
solver routine.  Whether the solver makes use of the solution in any
way is solver-dependent.
*/
void
OsiTMINLPInterface::setColSolution(const double *colsol)
{
  problem_->setxInit(getNumCols(), colsol);
  hasBeenOptimized_ = false;
  //  problem_->SetStartingPoint(getNumCols(), colsol);
}

/** Set dual solution variable values

rowprice[getNumRows()] is an array of values for the dual
variables. These values are copied to memory owned by the solver
interface object or the solver.  They will be returned as the result of
getRowPrice() until changed by another call to setRowPrice() or by a
call to any solver routine.  Whether the solver makes use of the
solution in any way is solver-dependent.
*/

void
OsiTMINLPInterface::setRowPrice(const double * rowprice)
{
  problem_->setDualsInit(getNumCols()*2 + getNumRows(), rowprice);
  hasBeenOptimized_ = false;
}

  /*! \brief Get an empty warm start object

  This routine returns an empty CoinWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
  */
CoinWarmStart *
OsiTMINLPInterface::getEmptyWarmStart () const
{return app_->getEmptyWarmStart();}

  /** Get warmstarting information */
CoinWarmStart* 
OsiTMINLPInterface::getWarmStart() const
{
  if (expose_warm_start_) {
    return app_->getWarmStart(problem_);
  }
  else {
    return getEmptyWarmStart();
  }
}
  /** Set warmstarting information. Return true/false depending on whether
      the warmstart information was accepted or not. */
bool 
OsiTMINLPInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  hasBeenOptimized_ = false;
  if (expose_warm_start_) {
    return app_->setWarmStart(warmstart, problem_);
  }
  else {
    return true;
  }
}

/** Set the index-th variable to be a continuous variable */
void
OsiTMINLPInterface::setContinuous(int index)
{
  problem_->SetVariableType(index, Bonmin::TMINLP::CONTINUOUS);
  hasBeenOptimized_ = false;
}
/** Set the index-th variable to be an integer variable */
void
OsiTMINLPInterface::setInteger(int index)
{
  problem_->SetVariableType(index, Bonmin::TMINLP::INTEGER);
  hasBeenOptimized_ = false;
}

/// Get objective function value (can't use default)
double
OsiTMINLPInterface::getObjValue() const
{
  return problem_->obj_value();
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiTMINLPInterface::setIntParam(OsiIntParam key, int value)
{
  //  debugMessage("OsiCpxSolverInterface::setIntParam(%d, %d)\n", key, value);

  bool retval = false;
  switch (key) {
  case OsiMaxNumIteration:
    retval = false;
    break;
  case OsiMaxNumIterationHotStart:
    if( value >= 0 ) {
      retval = false;
    }
    else
      retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiTMINLPInterface::setDblParam(OsiDblParam key, double value)
{
  //  debugMessage("OsiTMINLPInterface::setDblParam(%d, %g)\n", key, value);

  bool retval = false;
  switch (key) {
  case OsiDualObjectiveLimit:
    OsiDualObjectiveLimit_ = value;
    retval = true;
    break;
  case OsiPrimalObjectiveLimit:
    std::cerr<<"Can not set primal objective limit parameter"<<std::endl;
    retval = false;
    break;
  case OsiDualTolerance:
    std::cerr<<"Can not set dual tolerance parameter"<<std::endl;
    retval = false;
    break;
  case OsiPrimalTolerance:
    std::cerr<<"Can not set primal tolerance parameter"<<std::endl;
    retval = false;
  case OsiObjOffset:
    retval = OsiSolverInterface::setDblParam(key,value);
    break;
  case OsiLastDblParam:
    retval = false;
    break;
  }
  return retval;
}


//-----------------------------------------------------------------------------

bool
OsiTMINLPInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  //  debugMessage("OsiTMINLPInterface::setStrParam(%d, %s)\n", key, value.c_str());

  bool retval=false;
  switch (key) {
  case OsiProbName:
    OsiSolverInterface::setStrParam(key,value);
    return retval = true;
  case OsiSolverName:
    return false;
  case OsiLastStrParam:
    return false;
  }
  return false;
}

//-----------------------------------------------------------------------------

bool
OsiTMINLPInterface::getIntParam(OsiIntParam key, int& value) const
{
  //  debugMessage("OsiTMINLPInterface::getIntParam(%d)\n", key);

  bool retval = false;
  switch (key) {
  case OsiMaxNumIteration:
    retval = false;
    break;
  case OsiMaxNumIterationHotStart:
    retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiTMINLPInterface::getDblParam(OsiDblParam key, double& value) const
{
  //  debugMessage("OsiTMINLPInterface::getDblParam(%d)\n", key);

  bool retval = false;
  switch (key) {
  case OsiDualObjectiveLimit:
    value = OsiDualObjectiveLimit_;
    retval = true;
    break;
  case OsiPrimalObjectiveLimit:
    value = getInfinity();
    retval = true;
    break;
  case OsiDualTolerance:
    retval = false;
    break;
  case OsiPrimalTolerance:
    retval = false;
    break;
  case OsiObjOffset:
    retval = OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiLastDblParam:
    retval = false;
    break;
  }
  return retval;
}


//-----------------------------------------------------------------------------

bool
OsiTMINLPInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  //  debugMessage("OsiTMINLPInterface::getStrParam(%d)\n", key);

  switch (key) {
  case OsiProbName:
    OsiSolverInterface::getStrParam(key, value);
    break;
  case OsiSolverName:
    value = "Ipopt";
    break;
  case OsiLastStrParam:
    return false;
  }

  return true;
}

void
OsiTMINLPInterface::randomStartingPoint()
{
  int numcols = getNumCols();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  double * sol = new double[numcols];
  for(int i = 0 ; i < numcols ; i++) {
    double lower = min(-maxRandomRadius_,colUpper[i] - maxRandomRadius_);
    lower = max(colLower[i], lower);
    double upper = max(maxRandomRadius_,colLower[i] + maxRandomRadius_);
    upper = min(colUpper[i],upper);
    lower = min(upper,lower);
    upper = max(upper, lower);
    double interval = upper - lower;
    sol[i] = CoinDrand48()*(interval) + lower;
//    printf("%f in [%f,%f]\n",sol[i],lower,upper);
    //  std::cout<<interval<<"\t";
  }
  //std::cout<<std::endl;
  app_->disableWarmStart();
  setColSolution(sol);
  delete [] sol;
}



/** This methods initialiaze arrays for storing the jacobian */
int OsiTMINLPInterface::initializeJacobianArrays()
{
  Ipopt::Index n, m, nnz_h_lag;
  Ipopt::TNLP::IndexStyleEnum index_style;
  tminlp_->get_nlp_info( n, m, nnz_jac, nnz_h_lag, index_style);

  if(jRow_ != NULL) delete jRow_;
  if(jCol_ != NULL) delete jCol_;
  if(jValues_ != NULL) delete jValues_;

  jRow_ = new Ipopt::Index[nnz_jac];
  jCol_ = new Ipopt::Index[nnz_jac];
  jValues_ = new Ipopt::Number[nnz_jac];
  tminlp_->eval_jac_g(n, NULL, 0, m, nnz_jac, jRow_, jCol_, NULL);


  if(constTypes_ != NULL) delete [] constTypes_;
//  if(constTypesNum_ != NULL) delete [] constTypesNum_;

  constTypes_ = new Ipopt::TNLP::LinearityType[getNumRows()];
  tminlp_->get_constraints_linearity(getNumRows(), constTypes_);
//  constTypesNum_ = new int[getNumRows()];
  for(int i = 0; i < getNumRows() ; i++) {
    if(constTypes_[i]==Ipopt::TNLP::LINEAR) {
      //constTypesNum_[i] =
      nLinear_++;
    }
    else if(constTypes_[i]==Ipopt::TNLP::NON_LINEAR) {
      //constTypesNum_[i] = 
      nNonLinear_++;
    }
  }
  return nnz_jac;
}


/** get NLP constraint violation of point */
double
OsiTMINLPInterface::getConstraintViolation(const double *x, const double obj)
{
  int numcols = getNumCols();
  int numrows = getNumRows();
  double * g = new double[numrows];
  tminlp_->eval_g(numcols, x, 1, numrows, g);
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();


  double norm = 0;
  for(int i = 0; i< numrows ; i++) {
    if(constTypes_[i] == Ipopt::TNLP::NON_LINEAR) {
      double rowViolation = 0;
      if(rowLower[i] > -1e10)
         rowViolation = max(0.,rowLower[i] - g[i]);
//      if(rowViolation  > 0.) rowViolation  /= fabs(rowLower[i]);
      if(rowUpper[i] < 1e10);
        rowViolation = max(rowViolation, g[i] - rowUpper[i]);
//      if(rowViolation2 > 0.) rowViolation2 /= fabs(rowUpper[i]);
      norm = rowViolation > norm ? rowViolation : norm;
    }
  }

  double f;
  tminlp_->eval_f(numcols, x, 1, f);
  norm =  fabs(f - obj) > norm ? fabs(f - obj) : norm;

  delete [] g;
  return norm;
}


//A procedure to try to remove small coefficients in OA cuts (or make it non small
static inline
bool cleanNnz(double &value, double colLower, double colUpper,
    double rowLower, double rowUpper, double colsol,
    double & lb, double &ub, double tiny, double veryTiny)
{
  if(fabs(value)>= tiny) return 1;

  if(fabs(value)<veryTiny) return 0;//Take the risk?

  //try and remove
  double infty = 1e20;
  bool colUpBounded = colUpper < 10000;
  bool colLoBounded = colLower > -10000;
  bool rowNotLoBounded =  rowLower <= - infty;
  bool rowNotUpBounded = rowUpper >= infty;
  bool pos =  value > 0;

  if(!rowNotLoBounded && ! rowNotUpBounded)//would have to either choose side or duplicate cut
  {
    std::cerr<<"OA on non-convex constraint is very experimental, don't know how to remove"
    <<std::endl;
  }

  if(colLoBounded && pos && rowNotUpBounded) {
    lb += value * (colsol - colLower);
    return 0;
  }
  else
    if(colLoBounded && !pos && rowNotLoBounded) {
      ub += value * (colsol - colLower);
      return 0;
    }
    else
      if(colUpBounded && !pos && rowNotUpBounded) {
        lb += value * (colsol - colUpper);
        return 0;
      }
      else
        if(colUpBounded && pos && rowNotLoBounded) {
          ub += value * (colsol - colUpper);
          return 0;
        }
  //can not remove coefficient increase it to smallest non zero
  if(pos) value = tiny;
  else
    value = - tiny;
  return 1;
}

/** Get the outer approximation constraints at the current optimum of the
  current ipopt problem */
void
OsiTMINLPInterface::getOuterApproximation(OsiCuts &cs, bool getObj)
{
  getOuterApproximation(cs, getColSolution(), getObj);
}

/** Get the outer approximation constraints at the current optimum of the
  current ipopt problem */
void
OsiTMINLPInterface::getOuterApproximation(OsiCuts &cs, const double * x, bool getObj)
{
  int n,m, nnz_jac_g, nnz_h_lag;
  Ipopt::TNLP::IndexStyleEnum index_style;
  tminlp_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
  if(jRow_ == NULL || jCol_ == NULL || jValues_ == NULL)
    initializeJacobianArrays();
  assert(jRow_ != NULL);
  assert(jCol_ != NULL);
  double * g = new double[m];
  tminlp_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, jValues_);
  tminlp_->eval_g(n,x,1,m,g);
  //As jacobian is stored by cols fill OsiCuts with cuts
  CoinPackedVector * cuts = new CoinPackedVector[nNonLinear_ + 1];
  double * lb = new double[nNonLinear_ + 1];
  double * ub = new double[nNonLinear_ + 1];
  int * binding = new int[nNonLinear_ + 1];//store binding constraints at current opt (which are added to OA) -1 if constraint is not binding, otherwise index in subset of binding constraints
  int numBindings = 0;
    
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice();
  double infty = getInfinity();
  double nlp_infty = infty_;
  
  for(int i = 0; i< m ; i++) {
    if(constTypes_[i] == Ipopt::TNLP::NON_LINEAR) {
      if(rowLower[i] > - nlp_infty && rowUpper[i] < nlp_infty && fabs(duals[i]) == 0.)
      {
        binding[i] = -1;
        std::cerr<<"non binding constraint"<<std::endl;
        continue;
      }
      binding[i] = numBindings;
      if(rowLower[i] > - nlp_infty)
        lb[numBindings] = rowLower[i] - g[i];
      else
        lb[numBindings] = - infty;
      if(rowUpper[i] < nlp_infty)
        ub[numBindings] = rowUpper[i] - g[i];
      else
        ub[numBindings] = infty;
      if(rowLower[i] > -infty && rowUpper[i] < infty)
      {
        if(duals[i] >= 0)// <= inequality
          lb[numBindings] = - infty;
        if(duals[i] <= 0)// >= inequality
          ub[numBindings] = infty;
      }
      
      numBindings++;
    }
  }

  for(int i = 0 ; i < nnz_jac_g ; i++) {
    if(constTypes_[jRow_[i] - 1] == Ipopt::TNLP::NON_LINEAR) {
      //"clean" coefficient
      if(cleanNnz(jValues_[i],colLower[jCol_[i] - 1], colUpper[jCol_[i]-1],
          rowLower[jRow_[i] - 1], rowUpper[jRow_[i] - 1],
          x[jCol_[i] - 1],
          lb[binding[jRow_[i] - 1]],
          ub[binding[jRow_[i] - 1]], tiny_, veryTiny_)) {
        cuts[binding[jRow_[i] - 1]].insert(jCol_[i]-1,jValues_[i]);
        if(lb[binding[jRow_[i] - 1]] > - infty)
          lb[binding[jRow_[i] - 1]] += jValues_[i] * x[jCol_ [i] - 1];
        if(ub[binding[jRow_[i] - 1]] < infty)
        ub[binding[jRow_[i] - 1]] += jValues_[i] * x[jCol_ [i] - 1];
      }
    }
  }

  for(int i = 0; i< numBindings ; i++) {
    OsiRowCut newCut;
    //    if(lb[i]>-1e20) assert (ub[i]>1e20);

    newCut.setGloballyValid();
    newCut.setEffectiveness(99.99e99);
    newCut.setLb(lb[i]);
    newCut.setUb(ub[i]);
    newCut.setRow(cuts[i]);
    //    CftValidator validator;
    //    validator(newCut);
    cs.insert(newCut);
  }

  delete[] g;
  delete [] cuts;
  delete [] binding;

  if(getObj)  // Get the objective cuts
  {
    double * obj = new double [n];
    tminlp_->eval_grad_f(n, x, 1,obj);
    double f;
    tminlp_->eval_f(n, x, 1, f);

    CoinPackedVector v;
    v.reserve(n);
    lb[nNonLinear_] = -f;
    ub[nNonLinear_] = -f;
    //double minCoeff = 1e50;
    for(int i = 0; i<n ; i++)
    {
      if(cleanNnz(obj[i],colLower[i], colUpper[i],
          -getInfinity(), 0,
          x[i],
          lb[nNonLinear_],
          ub[nNonLinear_],tiny_, 1e-15)) {
        //	      minCoeff = min(fabs(obj[i]), minCoeff);
        v.insert(i,obj[i]);
        lb[nNonLinear_] += obj[i] * x[i];
        ub[nNonLinear_] += obj[i] * x[i];
      }
    }
    v.insert(n,-1);
    OsiRowCut newCut;
    newCut.setGloballyValid();
    newCut.setEffectiveness(99.99e99);
    newCut.setRow(v);
    newCut.setLb(-DBL_MAX/*Infinity*/);
    newCut.setUb(ub[nNonLinear_]);
//     CftValidator validator;
//     validator(newCut);
    cs.insert(newCut);
    delete [] obj;
  }

  delete []lb;
  delete[]ub;
}

double
OsiTMINLPInterface::getFeasibilityOuterApproximation(int n,const double * x_bar,const int *inds, OsiCuts &cs)
{
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_dist2point_obj(n,(const Ipopt::Number *) x_bar,(const Ipopt::Index *) inds);
  nCallOptimizeTNLP_++;
  totalNlpSolveTime_-=CoinCpuTime();
  Ipopt::SmartPtr<TNLPSolver> app2 = app_->clone();
  app2->Options()->SetIntegerValue("print_level", (Ipopt::Index) 0);
  optimization_status_ = app2->OptimizeTNLP(GetRawPtr(feasibilityProblem_));
  totalNlpSolveTime_+=CoinCpuTime();
  getOuterApproximation(cs, getColSolution(), 0);
  hasBeenOptimized_=true;
  return getObjValue();
}





void
OsiTMINLPInterface::extractLinearRelaxation(OsiSolverInterface &si, bool getObj,
                                            bool solveNlp)
{
  if(solveNlp)
    initialSolve();

  CoinPackedMatrix mat;
  double * rowLow = NULL;
  double * rowUp = NULL;

  int n;
  int m;
  int nnz_jac_g;
  int nnz_h_lag;
  Ipopt::TNLP::IndexStyleEnum index_style;
  //Get problem information
  tminlp_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);

  //if not allocated allocate spaced for stroring jacobian
  if(jRow_ == NULL || jCol_ == NULL || jValues_ == NULL)
    initializeJacobianArrays();

  //get Jacobian
  tminlp_->eval_jac_g(n, getColSolution(), 1, m, nnz_jac_g, NULL, NULL, jValues_);


  double *g = new double[m];
  tminlp_->eval_g(n, getColSolution(),1,m,g);

  rowLow = new double[m];
  rowUp = new double[m];
  int * binding = new int[m];//store binding constraints at current opt (which are added to OA) -1 if constraint is not binding, otherwise index in subset of binding constraints
  int numBindings = 0;
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice();
  assert(m==getNumRows());
  double infty = getInfinity();
  double nlp_infty = infty_;
  for(int i = 0 ; i < m ; i++) {
    {
      if(constTypes_[i] == Ipopt::TNLP::NON_LINEAR) {
        //If constraint is equality not binding do not add
        if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty && fabs(duals[i]) == 0.)
        {
            binding[i] = -1;
            continue;
        }
        else
          binding[i] = numBindings;
        if(rowLower[i] > - nlp_infty)
          rowLow[numBindings] = (rowLower[i] - g[i]) - 1e-07;
        else
          rowLow[numBindings] = - infty;
        if(rowUpper[i] < nlp_infty)
          rowUp[numBindings] =  (rowUpper[i] - g[i]) + 1e-07;
        else
          rowUp[numBindings] = infty;
        
        //If equality or ranged constraint only add one side by looking at sign of dual multiplier
        if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty)
        {
          if(duals[i] >= 0)// <= inequality
            rowLow[numBindings] = - infty;
          if(duals[i] <= 0)// >= inequality
            rowUp[numBindings] = infty;
        }
        numBindings++;
      }
      else {
        binding[i] = numBindings;//All are considered as bindings
        rowLow[numBindings] = (rowLower[i] - g[i]);
        rowUp[numBindings] =  (rowUpper[i] - g[i]);
        numBindings++;
      }
    }
  }

  //Then convert everything to a CoinPackedMatrix (col ordered)
  CoinBigIndex * inds = new CoinBigIndex[nnz_jac_g + 1];
  double * vals = new double [nnz_jac_g + 1];
  int * start = new int[n+1];
  int * length = new int[n];
  bool needOrder = false;
  int nnz = 0;
  if(nnz_jac_g > 0)
  {
  for(int k = 0; k < jCol_[0]; k++) {
    start[k] = nnz;
    length[k] = 0;
  }
  for(int i = 0 ; i < nnz_jac_g - 1 ; i++) {
    {
      if(jCol_[i + 1] < jCol_ [i] || ( jCol_ [i+1] == jCol_[i] && jRow_[i + 1] <= jRow_[i]) ) {
        needOrder = 1;
        break;
      }
      if(Ipopt::TNLP::LINEAR //Always accept coefficients from linear constraints
          || //For other clean tinys
          cleanNnz(jValues_[i],colLower[jCol_[i] - 1], colUpper[jCol_[i]-1],
              rowLower[jRow_[i] - 1], rowUpper[jRow_[i] - 1],
              getColSolution()[jCol_[i] - 1],
              rowLow[binding[jRow_[i] - 1]],
              rowUp[binding[jRow_[i] -1]], tiny_, veryTiny_)) {
        if(binding[jRow_[i] - 1] == -1) continue;
        vals[nnz] = jValues_[i];
        if(rowLow[binding[jRow_[i] - 1]] > - infty_)
          rowLow[binding[jRow_[i] - 1]] += jValues_[i] * getColSolution()[jCol_ [i] - 1];
        if(rowUp[binding[jRow_[i] - 1]] < infty_)
          rowUp[binding[jRow_[i] -1]] += jValues_[i] *getColSolution()[jCol_[i] -1];
        inds[nnz] = binding[jRow_[i] - 1];//jRow_[i ] - 1;
        length[jCol_[i] - 1]++;
        nnz++;
      }
    }
    if(jCol_[i + 1] > jCol_[i]) {
      for(int k = jCol_[i]; k < jCol_[i + 1]; k++) {
        start[k] = nnz;
        length[k] = 0;
      }
    }
  }
  if(!needOrder) {
    {
      length[jCol_[nnz_jac_g -1] - 1]++;
      vals[nnz] = jValues_[nnz_jac_g - 1];
      rowLow[binding[jRow_[nnz_jac_g - 1] - 1]] += jValues_[nnz_jac_g - 1] * getColSolution()[jCol_ [nnz_jac_g - 1] - 1];
      rowUp[binding[jRow_[nnz_jac_g - 1] -1]] += jValues_[nnz_jac_g - 1] *getColSolution()[jCol_[nnz_jac_g - 1] -1];
      inds[nnz++] = binding[jRow_[nnz_jac_g - 1] - 1];
    }
    for(int i = jCol_[nnz_jac_g -1] ; i < n ; i++) {
      start[i] = nnz;
      length[i] = 0;
    }
    start[n]=nnz;
  }
  else {
    std::cerr<<"jacobian matrix is not ordered properly"<<std::endl;
    throw -1;
  }
  delete [] g;
  }
  else {
   for (int i = 0 ; i < n ; i++)
   {
     length[i] = 0;
     start[i] = 0;
   }
   start[n]=0;
 }
 mat.assignMatrix(false, m, n, nnz, vals, inds, start, length);
  int numcols=getNumCols();
  double *obj = new double[numcols];
  for(int i = 0 ; i < numcols ; i++)
    obj[i] = 0;
  mat.transpose();
  
#if 0 
  std::cout<<"Mat is ordered by "<<
  <<mat.isColOrdered?"columns.":"rows."
  <<std::endl;
#endif
  
  si.loadProblem(mat, getColLower(), getColUpper(), obj, rowLow, rowUp);
  delete [] rowLow;
  delete [] rowUp;
  delete [] binding;
  for(int i = 0 ; i < getNumCols() ; i++) {
    if(isInteger(i))
      si.setInteger(i);
  }
  if(getObj) {
    //add variable alpha
    //(alpha should be empty in the matrix with a coefficient of -1 and unbounded)
    CoinPackedVector a;
    si.addCol(a,-si.getInfinity(), si.getInfinity(), 1.);

    // Now get the objective cuts
    // get the gradient, pack it and add the cut
    tminlp_->eval_grad_f(n, getColSolution(), 1,obj);
    double ub;
    tminlp_->eval_f(n, getColSolution(), 1, ub);
    ub*=-1;
    double lb = -1e300;
    CoinPackedVector objCut;
    CoinPackedVector * v = &objCut;
    v->reserve(n);
    for(int i = 0; i<n ; i++) {
     if(nnz_jac_g)
     {
      if(cleanNnz(obj[i],colLower[i], colUpper[i],
          -getInfinity(), 0,
          getColSolution()[i],
          lb,
          ub, tiny_, veryTiny_)) {
        v->insert(i,obj[i]);
        lb += obj[i] * getColSolution()[i];
        ub += obj[i] * getColSolution()[i];
      }
     }
     else //Unconstrained problem can not put clean coefficient
     {
         if(cleanNnz(obj[i],colLower[i], colUpper[i],
          -getInfinity(), 0,
          getColSolution()[i],
          lb,
          ub, 1e-03, 1e-08)) {
        v->insert(i,obj[i]);
        lb += obj[i] * getColSolution()[i];
        ub += obj[i] * getColSolution()[i];
         }
     }
    }
    v->insert(n,-1);
    si.addRow(objCut, lb, ub);
  }
  delete [] obj;
  
  if(solveNlp){
   app_->enableWarmStart();
   setColSolution(problem()->x_sol());
   setRowPrice(problem()->duals_sol());}

}

/** Add a collection of linear cuts to problem formulation.*/
void 
OsiTMINLPInterface::applyRowCuts(int numberCuts, const OsiRowCut * cuts)
{
  const OsiRowCut ** cutsPtrs = new const OsiRowCut*[numberCuts];
  for(int i = 0 ; i < numberCuts ; i++)
  {
    cutsPtrs[i] = &cuts[i];
  }
  tminlp_->addCuts(numberCuts, cutsPtrs);
  delete [] cutsPtrs;
}

void
OsiTMINLPInterface::solveAndCheckErrors(bool warmStarted, bool throwOnFailure,
    const char * whereFrom)
{
  totalNlpSolveTime_-=CoinCpuTime();
  if(warmStarted)
    optimization_status_ = app_->ReOptimizeTNLP(GetRawPtr(problem_));
  else
    optimization_status_ = app_->OptimizeTNLP(GetRawPtr(problem_));
  totalNlpSolveTime_+=CoinCpuTime();
  nCallOptimizeTNLP_++;
  hasBeenOptimized_ = true;

  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    //app_->Options()->SetIntegerValue("print_level",0, true, true);
    app_->Options()->SetStringValue("print_user_options","no", false, true);
  }


#if 1
  if(optimization_status_ == TNLPSolver::notEnoughFreedom)//Too few degrees of freedom
    {
      std::cout<<"Too few degrees of freedom...."<<std::endl;
      int numrows = getNumRows();
      int numcols = getNumCols();

      const double * colLower = getColLower();
      const double * colUpper = getColUpper();


      const double * rowLower = getRowLower();
      const double * rowUpper = getRowUpper();

      int numberFixed = 0;
      for(int i = 0 ; i < numcols ; i++)
	{
	  if(colUpper[i] - colLower[i] <= INT_BIAS)
	    {
	      numberFixed++;
	    }
	}
      int numberEqualities = 0;
      for(int i = 0 ; i < numrows ; i++)
	{
	  if(rowUpper[i] - rowLower[i] <= INT_BIAS)
	    {
	      numberEqualities++;
	    }	  
	}
      if(numcols - numberFixed > numberEqualities)
	{
	  std::string probName;
	  getStrParam(OsiProbName, probName);
	  throw newUnsolvedError(optimization_status_, problem_, probName);
	}
      double * saveColLow = CoinCopyOfArray(getColLower(), getNumCols());
      double * saveColUp = CoinCopyOfArray(getColUpper(), getNumCols());

      for(int i = 0; i < numcols && numcols - numberFixed <= numberEqualities ; i++)
	{
	  if(colUpper[i] - colLower[i] <= INT_BIAS)
	    {
	      setColLower(i, saveColLow[i]-1e-06);
	      setColUpper(i, saveColLow[i]+1e-06);
	      numberFixed--;
	    }
	}
      solveAndCheckErrors(warmStarted, throwOnFailure, whereFrom);
      //restore
      for(int i = 0; i < numcols && numcols - numberFixed < numrows ; i++)
	{
	  problem_->SetVariableLowerBound(i,saveColLow[i]);
	  problem_->SetVariableUpperBound(i,saveColUp[i]);
	}
      delete [] saveColLow;
      delete [] saveColUp;
      return;
    }
  else 
#endif
#ifdef AskPierreIfThisGoesOut
    if(optimization_status_ < -1)//Solver failed and the error can not be recovered, throw it
  {
    std::string probName;
    getStrParam(OsiProbName, probName);
    throw newUnsolvedError(optimization_status_, problem_, probName);
  }
#endif
  try{
  totalIterations_ += app_->IterationCount();
  }
  catch(SimpleError &E)
    {
      if (throwOnFailure)//something failed throw
	{
	  throw SimpleError("No statistics available from Ipopt",whereFrom);
	}
      else return;
    }

  messageHandler()->message(IPOPT_SUMMARY, messages_)
  <<whereFrom<<optimization_status_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;

  if((nCallOptimizeTNLP_ % 20) == 1)
    messageHandler()->message(LOG_HEAD, messages_)<<CoinMessageEol;


  if (numIterationSuspect_ >= 0 && (getIterationCount()>numIterationSuspect_ || isAbandoned())) {
    messageHandler()->message(SUSPECT_PROBLEM,
        messages_)<<nCallOptimizeTNLP_<<CoinMessageEol;
    std::string subProbName;
    getStrParam(OsiProbName, subProbName);
    std::ostringstream os;
    os<<"_"<<nCallOptimizeTNLP_;
    subProbName+=os.str();
    problem_->outputDiffs(subProbName, NULL/*getVarNames()*/);
  }

}

////////////////////////////////////////////////////////////////////
// Solve Methods                                                  //
////////////////////////////////////////////////////////////////////
/// Solve initial continuous relaxation
void OsiTMINLPInterface::initialSolve()
{
  assert(IsValid(app_));
  assert(IsValid(problem_));

  if(!hasPrintedOptions) {
    int printOptions;
    app_->Options()->GetEnumValue("print_user_options",printOptions,"bonmin.");
    if(printOptions)
      app_->Options()->SetStringValue("print_user_options","yes");
  }
  solveAndCheckErrors(0,1,"initialSolve");
  
  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    app_->Options()->SetStringValue("print_user_options","no");
    app_->Options()->SetIntegerValue("print_level",0);
  }
  
  const char * status=OPT_SYMB;
  ;
  if(isAbandoned()) status=FAILED_SYMB;
  else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;
  messageHandler()->message(LOG_FIRST_LINE, messages_)<<nCallOptimizeTNLP_
						      <<status<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;
  
  if(isAbandoned()) {
    resolveForRobustness(numRetryUnsolved_);
  }
  else if(numRetryInitial_)
    {
      resolveForCost(numRetryInitial_);
      /** Only do it once at the root.*/
      numRetryInitial_ = 0;
    }
  
}

/** Resolve the continuous relaxation after problem modification.
 * \note for Ipopt, same as resolve */
void
OsiTMINLPInterface::resolve()
{
  assert(IsValid(app_));
  assert(IsValid(problem_));
  if (INT_BIAS > 0.) {
    app_->Options()->SetStringValue("warm_start_same_structure", "yes");
  }
  else {
    app_->Options()->SetStringValue("warm_start_same_structure", "no");
  }

  app_->enableWarmStart();

  solveAndCheckErrors(1,1,"resolve");
  
  const char * status=OPT_SYMB;
  ;
  if(isAbandoned()) status=FAILED_SYMB;
  else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;
  messageHandler()->message(LOG_FIRST_LINE, messages_)<<nCallOptimizeTNLP_
						      <<status<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;
  
  if(isAbandoned()) {
    resolveForRobustness(numRetryUnsolved_);
  }
  else if(numRetryResolve_)
    resolveForCost(numRetryResolve_);
  
}


////////////////////////////////////////////////////////////////
// Methods returning info on how the solution process terminated  //
////////////////////////////////////////////////////////////////
/// Are there a numerical difficulties?
bool OsiTMINLPInterface::isAbandoned() const
{
  return (
        (optimization_status_==TNLPSolver::iterationLimit)||
        (optimization_status_==TNLPSolver::computationError)||
        (optimization_status_==TNLPSolver::illDefinedProblem)||
        (optimization_status_==TNLPSolver::illegalOption)||
        (optimization_status_==TNLPSolver::externalException)||
        (optimization_status_==TNLPSolver::exception)
      );
}

/// Is optimality proven?
bool OsiTMINLPInterface::isProvenOptimal() const
{
  return (optimization_status_==TNLPSolver::solvedOptimal ||
	  optimization_status_==TNLPSolver::solvedOptimalTol);
}
/// Is primal infeasiblity proven?
bool OsiTMINLPInterface::isProvenPrimalInfeasible() const
{
  return (optimization_status_ == TNLPSolver::provenInfeasible);
}
/// Is dual infeasiblity proven?
bool OsiTMINLPInterface::isProvenDualInfeasible() const
{
  throw SimpleError("Don't have this optimization status yet.",
      "isProvenDualInfeasible");
}
/// Is the given primal objective limit reached?
bool OsiTMINLPInterface::isPrimalObjectiveLimitReached() const
{
  std::cerr<<"Warning : isPrimalObjectiveLimitReached not implemented yet"<<std::endl;
  return 0;
}
/// Is the given dual objective limit reached?
bool OsiTMINLPInterface::isDualObjectiveLimitReached() const
{
  //  std::cerr<<"Warning : isDualObjectiveLimitReached not implemented yet"<<std::endl;
  return (optimization_status_==TNLPSolver::unbounded);

}
/// Iteration limit reached?
bool OsiTMINLPInterface::isIterationLimitReached() const
{
  return (optimization_status_==TNLPSolver::iterationLimit);
}

void
OsiTMINLPInterface::extractInterfaceParams()
{
  if (IsValid(app_)) {
#ifdef COIN_HAS_FILTERSQP
    FilterSolver * filter = dynamic_cast<FilterSolver *>(Ipopt::GetRawPtr(app_));
#endif
    app_->Options()->GetNumericValue("max_random_point_radius",maxRandomRadius_,"bonmin.");

#ifdef COIN_HAS_FILTERSQP
    if(filter && maxRandomRadius_ > 10.){
      std::cerr<<"filterSqp picking a random point with a big value seems not to work well setting max_random_point_radius to 10"<<std::endl;
      maxRandomRadius_ = 10.;
    }
#endif
    app_->Options()->GetIntegerValue("num_retry_unsolved_random_point", numRetryUnsolved_,"bonmin.");
    app_->Options()->GetIntegerValue("num_resolve_at_root", numRetryInitial_,"bonmin.");
    app_->Options()->GetIntegerValue("num_resolve_at_node", numRetryResolve_,"bonmin.");
    app_->Options()->GetIntegerValue("num_iterations_suspect", numIterationSuspect_,"bonmin.");
    app_->Options()->GetEnumValue("nlp_failure_behavior",pretendFailIsInfeasible_,"bonmin.");
    app_->Options()->GetNumericValue
    ("warm_start_bound_frac" ,pushValue_,"bonmin.");
    app_->Options()->GetNumericValue("tiny_element",tiny_,"bonmin.");
    app_->Options()->GetNumericValue("very_tiny_element",veryTiny_,"bonmin.");
  }
}


}/** end namespace Bonmin*/

