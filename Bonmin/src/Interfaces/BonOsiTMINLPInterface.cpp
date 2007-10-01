// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
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
#include "CoinTime.hpp"
#include <climits>
#include <string>
#include <sstream>
#include "BonAuxInfos.hpp"

#include "Ipopt/BonIpoptSolver.hpp"
#ifdef COIN_HAS_FILTERSQP
#include "Filter/BonFilterSolver.hpp"
#endif

#include "OsiBranchingObject.hpp"
#include "BonStrongBranchingSolver.hpp"

using namespace Ipopt;


namespace Bonmin {
///Register options
static void
register_general_options
(SmartPtr<RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("bonmin nlp interface option", RegisteredOptions::BonminCategory);
  roptions->AddStringOption2("nlp_solver",
                             "Choice of the solver for local optima of continuous nlp's",
                             "Ipopt",
                             "Ipopt", "Interior Point OPTimizer (https://projects.coin-or.org/Ipopt)",
                             "filterSQP", "Sequential quadratic programming trust region algorithm (http://www-unix.mcs.anl.gov/~leyffer/solvers.html)",
                             "");
  roptions->setOptionExtraInfo("nlp_solver",15);
  roptions->AddBoundedIntegerOption("nlp_log_level",
      "specify NLP solver interface log level (independent from ipopt print_level).",
      0,2,1,
      "Set the level of output of the OsiTMINLPInterface : "
      "0 - none, 1 - normal, 2 - verbose"
                                   );
  roptions->setOptionExtraInfo("nlp_log_level",15);

  roptions->AddStringOption2("file_solution",
       "Write a file bonmin.sol with the solution",
       "no",
       "yes","",
       "no","","");

  roptions->AddStringOption3("warm_start",
      "Select the warm start method",
      "none",
      "none","No warm start",
      "optimum","Warm start with direct parent optimum",
      "interior_point","Warm start with an interior point of direct parent",
      "This will affect the function getWarmStart(), and as a consequence the wam starting in the various algorithms.");
  roptions->setOptionExtraInfo("warm_start",8);

  roptions->SetRegisteringCategory("bonmin options for robustness", RegisteredOptions::BonminCategory);

  roptions->AddLowerBoundedNumberOption("max_random_point_radius",
      "Set max value r for coordinate of a random point.",
      0.,1,1e5,
      "When picking a random point coordinate i will be in the interval [min(max(l,-r),u-r), max(min(u,r),l+r)] "
      "(where l is the lower bound for the variable and u is its upper bound)");
  roptions->setOptionExtraInfo("max_random_point_radius",8);

  roptions->AddStringOption3("random_point_type","method to choose a random starting point",
			     "Jon",
			     "Jon", "Choose random point uniformly between the bounds",
			     "Andreas", "perturb the starting point of the problem within a prescriped interval",
			     "Claudia", "perturb the starting point using the perturbation radius suffix information",
			     "");
  roptions->setOptionExtraInfo("random_point_type",8);

    roptions->AddLowerBoundedNumberOption("random_point_perturbation_interval",
					   "Amount by which starting point is perturbed when choosing to pick random point by perturbating starting point",
					   0.,true, 1.,
					   "");
  roptions->setOptionExtraInfo("random_point_perturbation_interval",8);
					   

  roptions->AddLowerBoundedIntegerOption
  ("num_iterations_suspect",
   "Number of iterations over which a node is considered \"suspect\" (for debugging purposes only, see detailed documentation).",
   -1,-1,
   "When the number of iterations to solve a node is above this number, the subproblem at this"
   " node is considered to be suspect and it will be outputed in a file (set to -1 to deactivate this).");
  roptions->setOptionExtraInfo("num_iterations_suspect",15);

  

  roptions->AddLowerBoundedIntegerOption("num_retry_unsolved_random_point",
      "Number $k$ of times that the algorithm will try to resolve an unsolved NLP with a random starting point "
      "(we call unsolved an NLP for which Ipopt is not "
      "able to guarantee optimality within the specified tolerances).",
      0,0,
      "When Ipopt fails to solve a continuous NLP sub-problem, if $k > 0$, the algorithm will "
      "try again to solve the failed NLP with $k$ new randomly chosen starting points "
      " or until the problem is solved with success.");
  roptions->setOptionExtraInfo("num_retry_unsolved_random_point",15);


  roptions->SetRegisteringCategory("bonmin options for non-convex problems", RegisteredOptions::BonminCategory);


  roptions->AddLowerBoundedIntegerOption("num_resolve_at_root",
      "Number $k$ of tries to resolve the root node with different starting points.",
      0,0,
      "The algorithm will solve the root node with $k$ random starting points"
      " and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_root",8);

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_node",
      "Number $k$ of tries to resolve a node (other than the root) of the tree with different starting point.",
      0,0,
      "The algorithm will solve all the nodes with $k$ different random starting points "
      "and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_node",8);

  roptions->AddLowerBoundedIntegerOption("num_resolve_at_infeasibles",
      "Number $k$ of tries to resolve an infeasible node (other than the root) of the tree with different starting point.",
      0,0,
      "The algorithm will solve all the infeasible nodes with $k$ different random starting points "
      "and will keep the best local optimum found.");
  roptions->setOptionExtraInfo("num_resolve_at_infeasibles",8);



  }

static void register_OA_options
(SmartPtr<RegisteredOptions> roptions)
{
  roptions->SetRegisteringCategory("bonmin options : Outer Approximation cuts", RegisteredOptions::BonminCategory);
  
  roptions->AddStringOption2("disjunctive_cut_type",
      "Determine if and what kind of disjunctive cuts should be computed.",
      "none",
      "none", "No disjunctive cuts.",
      "most-fractional", "If discrete variables present, compute disjunction for most-fractional variable");
  roptions->setOptionExtraInfo("disjunctive_cut_type",7);

  roptions->AddStringOption2("oa_cuts_scope","Specify if OA cuts added are to be set globally or locally valid",
                             "global",
			     "local","Cuts are treated as globally valid",
			     "global", "Cuts are treated as locally valid",
			     "");
  roptions->setOptionExtraInfo("oa_cuts_scope",7);

  roptions->AddStringOption2("add_only_violated_oa","Do we add all OA cuts or only the ones violated by current point?",
			     "no",
			     "no","Add all cuts",
			     "yes","Add only violated Cuts","");
  roptions->setOptionExtraInfo("add_only_violated_oa",7);

  roptions->AddStringOption4("cut_strengthening_type",
                             "Determines if and what kind of cut strengthening should be performed.",
                             "none",
                             "none", "No strengthening of cuts.",
                             "sglobal", "Strengthen global cuts.",
                             "uglobal-slocal", "Unstrengthened global and strengthened local cuts",
                             "sglobal-slocal", "Strengthened global and strengthened local cuts",
                             "");
  roptions->setOptionExtraInfo("cut_strengthening_type",7);
  
  roptions->AddLowerBoundedNumberOption("tiny_element","Value for tiny element in OA cut",
      -0.,0,1e-08,
      "We will remove \"cleanly\" (by relaxing cut) an element lower"
      " than this.");
  roptions->setOptionExtraInfo("tiny_element",7);

  roptions->AddLowerBoundedNumberOption("very_tiny_element","Value for very tiny element in OA cut",
      -0.,0,1e-17,
      "Algorithm will take the risk of neglecting an element lower"
      " than this.");
  roptions->setOptionExtraInfo("very_tiny_element",7);

  roptions->AddLowerBoundedIntegerOption("oa_cuts_log_level",
                                         "level of log when generating OA cuts.",
                                         0, 0,
                                         "0: ouputs nothings,\n"
                                         "1: when a cut is generated, its violation and index of row from which it originates,\n"
                                         "2: always output violation of the cut.\n"
                                         "3: ouput generated cuts incidence vectors.");
  roptions->setOptionExtraInfo("oa_cuts_log_level",7);

}


///Register options
void
OsiTMINLPInterface::registerOptions
(SmartPtr<RegisteredOptions> roptions)
{
  // We try to register the options - if those have been registered
  // already, we catch the exception and don't need to do it again
  try {
    register_general_options(roptions);
    register_OA_options(roptions);
#ifdef COIN_HAS_FILTERSQP
    FilterSolver::RegisterOptions(roptions);
#endif
#ifdef COIN_HAS_IPOPT
    IpoptSolver::RegisterOptions(roptions);
#endif
  }   
  catch(RegisteredOptions::OPTION_ALREADY_REGISTERED) {
    // skipping
  }


}

/** To set some application specific defaults. */
void 
OsiTMINLPInterface::setAppDefaultOptions(SmartPtr<OptionsList> Options)
{}

OsiTMINLPInterface::Messages::Messages
():CoinMessages((int)OSITMINLPINTERFACE_DUMMY_END)
{
  strcpy(source_ ,"NLP");
  addMessage(SOLUTION_FOUND, CoinOneMessage
      (SOLUTION_FOUND + 1,2,"After %d tries found a solution of %g "
       "(previous best %g)."));

  addMessage(INFEASIBLE_SOLUTION_FOUND, CoinOneMessage
      (INFEASIBLE_SOLUTION_FOUND + 1 ,2,"After %d tries found an solution of %g "
       "infeasible problem."));

  addMessage(UNSOLVED_PROBLEM_FOUND, CoinOneMessage
      (UNSOLVED_PROBLEM_FOUND + 1,2,"After %d tries found an solution of %g "
       "unsolved problem."));
  addMessage(WARN_SUCCESS_WS,CoinOneMessage
      (WARN_SUCCESS_WS + 3001,2,
       "Problem not solved with warm start but "
       "solved without"));

  addMessage(WARNING_RESOLVING,CoinOneMessage
      (WARNING_RESOLVING + 3001,2,
       "Trying to resolve NLP with different starting "
       "point (%d attempts)."));
  addMessage(WARN_SUCCESS_RANDOM, CoinOneMessage
      (WARN_SUCCESS_RANDOM + 3001,1,
       "Problem initially not solved but solved with "
       "a random starting point (success on %d attempt)"));
  addMessage(WARN_CONTINUING_ON_FAILURE, CoinOneMessage
      (WARN_CONTINUING_ON_FAILURE + 3001, 1,
       "Warning : continuing branching, while there are "
       "unrecovered failures at nodes"));

  addMessage(SUSPECT_PROBLEM, CoinOneMessage
      (SUSPECT_PROBLEM + 1,2,"NLP number %d is suspect (see bounds and start file)"));
  addMessage(IPOPT_SUMMARY, CoinOneMessage
      (IPOPT_SUMMARY + 1,2,"Ipopt return (for %s): status %2d, iter count %4d, time %g"));
  addMessage(BETTER_SOL, CoinOneMessage
      (BETTER_SOL + 1,2,"Solution of value %g found on %d'th attempt"));

  addMessage(LOG_HEAD, CoinOneMessage
      (LOG_HEAD + 1,1,"\n          "
       "    Num      Status      Obj             It       time"));
  addMessage(LOG_FIRST_LINE, CoinOneMessage
      (LOG_FIRST_LINE + 1,1,
       "    %-8d %-11s %-14g %-8d %-3g"));
  addMessage(LOG_LINE, CoinOneMessage
      (LOG_LINE + 1,1,
       " %c  r%-7d %-11s %-14g %-8d %-3g"));

  addMessage(ALTERNATE_OBJECTIVE, CoinOneMessage
             (ALTERNATE_OBJECTIVE + 1,1,"Objective value recomputed with alternate objective: %g."));
  
  addMessage(WARN_RESOLVE_BEFORE_INITIAL_SOLVE, CoinOneMessage
      (WARN_RESOLVE_BEFORE_INITIAL_SOLVE + 3001,1,"resolve called before any call to initialSolve"
       " can not use warm starts."));
  addMessage(ERROR_NO_TNLPSOLVER, CoinOneMessage
             (ERROR_NO_TNLPSOLVER + 6001,1,"Can not parse options when no IpApplication has been created"));
  addMessage(WARNING_NON_CONVEX_OA, CoinOneMessage
             (WARNING_NON_CONVEX_OA + 3001, 1,
              "OA on non-convex constraint is very experimental."));                                          

}


void  
OsiTMINLPInterface::OaMessageHandler::print(OsiRowCut &row){
  FILE * fp = filePointer();
  const int &n = row.row().getNumElements();
  fprintf(fp,"Row cut has %d elements. Lower bound: %g, upper bound %g.\n", n, row.lb(), row.ub());
  const int * idx = row.row().getIndices();
  const double * val = row.row().getElements();
  for(int i = 0 ; i < n ; i++){
    fprintf(fp,"%g, x%d",val[i], idx[i]);
    if(i && i % 7 == 0)
      fprintf(fp,"\n");
  } 
}

OsiTMINLPInterface::OaMessages::OaMessages(): CoinMessages((int) OA_MESSAGES_DUMMY_END){
   strcpy(source_,"OaCg");
   addMessage(VIOLATED_OA_CUT_GENERATED, CoinOneMessage(
              VIOLATED_OA_CUT_GENERATED + 1, 1,"Row %d, cut violation is %g: Outer approximation cut generated."));

   addMessage(CUT_NOT_VIOLATED_ENOUGH, CoinOneMessage(
              CUT_NOT_VIOLATED_ENOUGH + 1, 2,"Row %d, cut violation is %g: Outer approximation cut not generated."));

   addMessage(OA_CUT_GENERATED, CoinOneMessage(
              OA_CUT_GENERATED + 1, 1,"Row %d: Outer approximation cut not generated."));
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
    hasVarNamesFile_(true),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(1e08),
    randomGenerationType_(0),
    max_perturbation_(COIN_DBL_MAX),
    pushValue_(1e-02),
    numRetryInitial_(-1),
    numRetryResolve_(-1),
    numRetryInfeasibles_(-1),
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
    nLinear_(0),
    nNonLinear_(0),
    tiny_(1e-8),
    veryTiny_(1e-20),
    infty_(1e100),
    exposeWarmStart_(false),
    firstSolve_(true),
    cutStrengthener_(NULL),
    oaMessages_(),
    oaHandler_(NULL)
{
   oaHandler_ = new OaMessageHandler;
   oaHandler_->setLogLevel(0);
}

void 
OsiTMINLPInterface::initialize(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                               Ipopt::SmartPtr<Ipopt::OptionsList> options,
                               Ipopt::SmartPtr<Ipopt::Journalist> journalist,
                               Ipopt::SmartPtr<TMINLP> tminlp){
  if(!IsValid(app_))
     createApplication(roptions, options, journalist);
  setModel(tminlp); 
}

void OsiTMINLPInterface::setSolver(Ipopt::SmartPtr<TNLPSolver> app){
  app_ = app;
  }


void
OsiTMINLPInterface::createApplication(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                                      Ipopt::SmartPtr<Ipopt::OptionsList> options,
                                      Ipopt::SmartPtr<Ipopt::Journalist> journalist
                                      )
{
  assert(!IsValid(app_));
  int ival;
  options->GetEnumValue("nlp_solver", ival, "bonmin.");
  Solver s = (Solver) ival;
  if(s == EFilterSQP){
#ifdef COIN_HAS_FILTERSQP
    app_ = new Bonmin::FilterSolver(roptions, options, journalist);
#else
#endif    
  }
  else if(s == EIpopt){
#ifdef COIN_HAS_IPOPT
    app_ = new IpoptSolver(roptions, options, journalist);
#else
    std::cerr<<"Bonmin not configured to run with Ipopt"<<std::endl;
    throw -1;
#endif
  }
  if (!app_->Initialize("")) {
    std::cerr<<"OsiTMINLPInterface: Error during initialization of app_"<<std::endl;
    throw -1;
  }
  extractInterfaceParams();
  
}

/** Facilitator to allocate a tminlp and an application. */
void
OsiTMINLPInterface::setModel(SmartPtr<TMINLP> tminlp)
{
  assert(IsValid(tminlp));
  tminlp_ = tminlp;
  problem_ = new TMINLP2TNLP(tminlp_);
}



void
OsiTMINLPInterface::readOptionFile(const std::string & fileName)
{
  if(IsValid(app_)){
  std::ifstream is;
  if (fileName != "") {
    try {
      is.open(fileName.c_str());
    }
    catch(std::bad_alloc) {
      std::cerr<<"Not enough memory to open option file.\n";
      throw -1;
    }
  }
  options()->ReadFromStream(*app_->journalist(), is);
  extractInterfaceParams();
  }
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
    hasVarNamesFile_(source.hasVarNamesFile_),
    nCallOptimizeTNLP_(0),
    totalNlpSolveTime_(0),
    totalIterations_(0),
    maxRandomRadius_(source.maxRandomRadius_),
    randomGenerationType_(source.randomGenerationType_),
    max_perturbation_(source.max_perturbation_),
    pushValue_(source.pushValue_),
    numRetryInitial_(source.numRetryInitial_),
    numRetryResolve_(source.numRetryResolve_),
    numRetryInfeasibles_(source.numRetryInfeasibles_),
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
    exposeWarmStart_(source.exposeWarmStart_),
    firstSolve_(true),
    cutStrengthener_(source.cutStrengthener_),
    oaMessages_(),
    oaHandler_(NULL),
    strong_branching_solver_(source.strong_branching_solver_)
{
      messageHandler()->setLogLevel(source.messageHandler()->logLevel());
  //Pass in message handler
  if(source.messageHandler())
    passInMessageHandler(source.messageHandler());
  // Copy options from old application
  if(IsValid(source.tminlp_)) {
    problem_ = new TMINLP2TNLP(tminlp_);
    problem_->copyUserModification(*source.problem_);
    pretendFailIsInfeasible_ = source.pretendFailIsInfeasible_;

    setAuxiliaryInfo(source.getAuxiliaryInfo());
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
    feasibilityProblem_ = new TNLP2FPNLP
        (SmartPtr<TNLP>(GetRawPtr(problem_)));
  else
    throw SimpleError("Don't know how to copy an empty IpoptOAInterface.",
        "copy constructor");


  if(source.jValues_!=NULL && source.jRow_ != NULL && source.jCol_ != NULL && nnz_jac>0) {
    jValues_ = new double [nnz_jac];
    jCol_    = new Index [nnz_jac];
    jRow_    = new Index [nnz_jac];
    CoinCopyN(source.jValues_ , nnz_jac,jValues_ );
    CoinCopyN(source.jCol_    , nnz_jac,jCol_    );
    CoinCopyN(source.jRow_    , nnz_jac,jRow_    );

    if(source.constTypes_ != NULL) {
      constTypes_ = new TNLP::LinearityType [getNumRows()];
      CoinCopyN(source.constTypes_, getNumRows(), constTypes_);
    }
  }
  else if(nnz_jac > 0) {
    throw SimpleError("Arrays for storing jacobian are inconsistant.",
        "copy constructor");
  }


   oaHandler_ = new OaMessageHandler(*source.oaHandler_);;

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
      problem_ = new TMINLP2TNLP(tminlp_);
      app_ = rhs.app_->clone();

      feasibilityProblem_ = new TNLP2FPNLP
          (SmartPtr<TNLP>(GetRawPtr(problem_)));
      nnz_jac = rhs.nnz_jac;

      if(constTypes_ != NULL) {
        delete [] constTypes_;
        constTypes_ = NULL;
      }
      if(rhs.constTypes_ != NULL) {
        constTypes_ = new TNLP::LinearityType[getNumRows()];
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
        jCol_    = new Index [nnz_jac];
        jRow_    = new Index [nnz_jac];
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
      exposeWarmStart_ = rhs.exposeWarmStart_;

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

    hasVarNamesFile_ = rhs.hasVarNamesFile_;

    nCallOptimizeTNLP_ = rhs.nCallOptimizeTNLP_;
    totalNlpSolveTime_ = rhs.totalNlpSolveTime_;
    totalIterations_ = rhs.totalIterations_;
    maxRandomRadius_ = rhs.maxRandomRadius_;
    pushValue_ = rhs.pushValue_;
    numRetryInitial_ = rhs.numRetryInitial_;
    numRetryResolve_ = rhs.numRetryResolve_;
    numRetryInfeasibles_ = rhs.numRetryInfeasibles_;
    numRetryUnsolved_ = rhs.numRetryUnsolved_;
    pretendFailIsInfeasible_ = rhs.pretendFailIsInfeasible_;
    numIterationSuspect_ = rhs.numIterationSuspect_;

    hasBeenOptimized_ = rhs.hasBeenOptimized_;
    cutStrengthener_ = rhs.cutStrengthener_;

    delete oaHandler_;
    oaHandler_ = new OaMessageHandler(*rhs.oaHandler_);
    strong_branching_solver_ = rhs.strong_branching_solver_;

    freeCachedData();
  }
  return *this;
}

SmartPtr<OptionsList> OsiTMINLPInterface::options()
{
  if(!IsValid(app_)) {
    messageHandler()->message(ERROR_NO_TNLPSOLVER, messages_)<<CoinMessageEol;
    return NULL;
  }
  else
    return app_->options();
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
  delete oaHandler_;
}

void
OsiTMINLPInterface::freeCachedColRim()
{
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
  setColSolution(point);
  setRowPrice(&point[getNumCols()]);
  app_->enableWarmStart();
  delete [] point;

  optimizationStatus_ = app_->ReOptimizeTNLP(GetRawPtr(problem_));
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
    <<"resolveForRobustness"<<optimizationStatus_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;


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
    throw newUnsolvedError(app_->errorCode(), problem_,
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

#if 1


///get name of variables
const OsiSolverInterface::OsiNameVec& 
OsiTMINLPInterface::getVarNames() {
  return getColNames();
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
  return (problem_->var_types()[colNumber]==TMINLP::CONTINUOUS);
}

/// Return true if column is binary
bool
OsiTMINLPInterface::isBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==TMINLP::BINARY);
}

/** Return true if column is integer.
    Note: This function returns true if the the column
    is binary or a general integer.
*/
bool
OsiTMINLPInterface::isInteger(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==TMINLP::BINARY)||
      (problem_->var_types()[colNumber]==TMINLP::INTEGER));
}

/// Return true if column is general integer
bool
OsiTMINLPInterface::isIntegerNonBinary(int colNumber) const
{
  return (problem_->var_types()[colNumber]==TMINLP::INTEGER);
}
/// Return true if column is binary and not fixed at either bound
bool
OsiTMINLPInterface::isFreeBinary(int colNumber) const
{
  return ((problem_->var_types()[colNumber]==TMINLP::BINARY)
      &&((getColUpper()[colNumber]-getColLower()[colNumber]) > 1 - 1e-09));
}

/// Get solver's value for infinity
double
OsiTMINLPInterface::getInfinity() const
{
  return COIN_DBL_MAX;
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
  if (exposeWarmStart_) {
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
  if (exposeWarmStart_) {
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
  problem_->SetVariableType(index, TMINLP::CONTINUOUS);
  hasBeenOptimized_ = false;
}
/** Set the index-th variable to be an integer variable */
void
OsiTMINLPInterface::setInteger(int index)
{
  problem_->SetVariableType(index, TMINLP::INTEGER);
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
  default:
    retval = false;
    std::cerr << "Unhandled case in setIntParam\n";
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
  default:
    retval = false;
    std::cerr << "Unhandled case in setDblParam\n";
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

  value = -COIN_INT_MAX; // Give a dummy value
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
  default:
    retval = false;
    std::cerr << "Unhandled case in setIntParam\n";
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
  const Number * x_init = problem_->x_init_user();
  const double* perturb_radius = NULL;
  if (randomGenerationType_ == perturb_suffix) {
    const TMINLP::PerturbInfo* pertubinfo = tminlp_->perturbInfo();
    if (pertubinfo) {
      perturb_radius = pertubinfo->GetPerturbationArray();
    }
    if (!perturb_radius) {
      throw SimpleError("Can't use perturb_radius if no radii are given.",
			"randomStartingPoint");
    }
  }
  for(int i = 0 ; i < numcols ; i++) {
    int randomGenerationType = randomGenerationType_;
    if(x_init[i] < colLower[i] || x_init[i] > colUpper[i])
      randomGenerationType = uniform;
    if(randomGenerationType_ == uniform){
      double lower = min(-maxRandomRadius_,colUpper[i] - maxRandomRadius_);
      lower = max(colLower[i], lower);
      double upper = max(maxRandomRadius_,colLower[i] + maxRandomRadius_);
      upper = min(colUpper[i],upper);
      lower = min(upper,lower);
      upper = max(upper, lower);
      double interval = upper - lower;
      sol[i] = CoinDrand48()*(interval) + lower;}
    else if (randomGenerationType_ == perturb){
      const double lower = max(x_init[i] - max_perturbation_, colLower[i]);
      const double upper = min(x_init[i] + max_perturbation_, colUpper[i]);
      const double interval = upper - lower;
      sol[i]  = lower + CoinDrand48()*(interval);
    }
    else if (randomGenerationType_ == perturb_suffix){
      const double radius = perturb_radius[i];
      const double lower = max(x_init[i] - radius*max_perturbation_, colLower[i]);
      const double upper = min(x_init[i] + radius*max_perturbation_, colUpper[i]);
      const double interval = upper - lower;
      sol[i]  = lower + CoinDrand48()*(interval);
    }
    //printf("%f in [%f,%f]\n",sol[i],lower,upper);
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
  Index n, m, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  tminlp_->get_nlp_info( n, m, nnz_jac, nnz_h_lag, index_style);

  if(jRow_ != NULL) delete jRow_;
  if(jCol_ != NULL) delete jCol_;
  if(jValues_ != NULL) delete jValues_;

  jRow_ = new Index[nnz_jac];
  jCol_ = new Index[nnz_jac];
  jValues_ = new Number[nnz_jac];
  tminlp_->eval_jac_g(n, NULL, 0, m, nnz_jac, jRow_, jCol_, NULL);
  if(index_style == Ipopt::TNLP::FORTRAN_STYLE)//put C-style
  {
    for(int i = 0 ; i < nnz_jac ; i++){
      jRow_[i]--;
      jCol_[i]--;
    }
  }

  if(constTypes_ != NULL) delete [] constTypes_;
//  if(constTypesNum_ != NULL) delete [] constTypesNum_;

  constTypes_ = new TNLP::LinearityType[getNumRows()];
  tminlp_->get_constraints_linearity(getNumRows(), constTypes_);
//  constTypesNum_ = new int[getNumRows()];
  for(int i = 0; i < getNumRows() ; i++) {
    if(constTypes_[i]==TNLP::LINEAR) {
      //constTypesNum_[i] =
      nLinear_++;
    }
    else if(constTypes_[i]==TNLP::NON_LINEAR) {
      //constTypesNum_[i] = 
      nNonLinear_++;
    }
  }
  return nnz_jac;
}


double 
OsiTMINLPInterface::getConstraintsViolation(const double *x, double &obj)
{
  int numcols = getNumCols();
  int numrows = getNumRows();
  double * g = new double[numrows];
  tminlp_->eval_g(numcols, x, 1, numrows, g);
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();


  double norm = 0;
  for(int i = 0; i< numrows ; i++) {
    if(!constTypes_ || constTypes_[i] == TNLP::NON_LINEAR) {
      double rowViolation = 0;
      if(rowLower[i] > -1e10)
         rowViolation = max(0.,rowLower[i] - g[i]);

      if(rowUpper[i] < 1e10);
        rowViolation = max(rowViolation, g[i] - rowUpper[i]);

      norm = rowViolation > norm ? rowViolation : norm;
    }
  }
  tminlp_->eval_f(numcols, x, 1, obj);
  delete [] g;
  return norm;
}

/** get infinity norm of constraint violation of a point and objective error*/
double
OsiTMINLPInterface::getNonLinearitiesViolation(const double *x, const double obj)
{
  double f;
  double norm = getConstraintsViolation(x, f);
  assert((f - obj) > -1e-08);
  norm =  (f - obj) > norm ? f - obj : norm;
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

/** Get the outer approximation constraints at the point x.
*/
void
OsiTMINLPInterface::getOuterApproximation(OsiCuts &cs, const double * x, bool getObj, const double * x2, double theta, bool global)
{
  int n,m, nnz_jac_g, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
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

  int * row2cutIdx = new int[m];//store correspondance between index of row and index of cut (some cuts are not generated for rows because linear, or not binding). -1 if constraint does not generate a cut, otherwise index in cuts.
  int numCuts = 0;

  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice() + 2 * n;
  double infty = getInfinity();
  double nlp_infty = infty_;
  
  for(int rowIdx = 0; rowIdx < m ; rowIdx++) {
    if(constTypes_[rowIdx] == TNLP::NON_LINEAR) {
#if 0
      if(fabs(duals[rowIdx]) == 0.)
      {
        row2cutIdx[rowIdx] = -1;
#ifdef NDEBUG
        std::cerr<<"non binding constraint"<<std::endl;
#endif
        continue;
      }
#endif
      row2cutIdx[rowIdx] = numCuts;
      if(rowLower[rowIdx] > - nlp_infty)
        lb[numCuts] = rowLower[rowIdx] - g[rowIdx];
      else
        lb[numCuts] = - infty;
      if(rowUpper[rowIdx] < nlp_infty)
        ub[numCuts] = rowUpper[rowIdx] - g[rowIdx];
      else
        ub[numCuts] = infty;
      if(rowLower[rowIdx] > -infty && rowUpper[rowIdx] < infty)
      {
        if(duals[rowIdx] >= 0)// <= inequality
          lb[numCuts] = - infty;
        if(duals[rowIdx] <= 0)// >= inequality
          ub[numCuts] = infty;
      }
      
      numCuts++;
    }
    else
      row2cutIdx[rowIdx] = -1;
  }


  for(int i = 0 ; i < nnz_jac_g ; i++) {
    const int &rowIdx = jRow_[i];
    const int & cutIdx = row2cutIdx[ rowIdx ];
    if(cutIdx != -1) {
      const int &colIdx = jCol_[i];
      //"clean" coefficient
      if(cleanNnz(jValues_[i],colLower[colIdx], colUpper[colIdx],
		  rowLower[rowIdx], rowUpper[rowIdx],
		  x[colIdx],
		  lb[cutIdx],
		  ub[cutIdx], tiny_, veryTiny_)) {
        cuts[cutIdx].insert(colIdx,jValues_[i]);
        if(lb[cutIdx] > - infty)
          lb[cutIdx] += jValues_[i] * x[colIdx];
        if(ub[cutIdx] < infty)
	  ub[cutIdx] += jValues_[i] * x[colIdx];
      }
    }
  }

  int * cut2rowIdx = NULL;
  if (IsValid(cutStrengthener_) || oaHandler_->logLevel() > 0) {
    cut2rowIdx = new int [numCuts];// Store correspondance between indices of cut and indices of rows. For each cut
    for(int rowIdx = 0 ; rowIdx < m ; rowIdx++){
       if(row2cutIdx[rowIdx] >= 0){
          cut2rowIdx[row2cutIdx[rowIdx]] = rowIdx;
       }
    }
  }

  for(int cutIdx = 0; cutIdx < numCuts; cutIdx++) {
    //Compute cut violation
    if(x2 != NULL) {
      double rhs = cuts[cutIdx].dotProduct(x2);
      double violation = 0.;
      violation = max(violation, rhs - ub[cutIdx]);
      violation = max(violation, lb[cutIdx] - rhs);
      if(violation < theta) {
        if(oaHandler_->logLevel() > 0)
          oaHandler_->message(CUT_NOT_VIOLATED_ENOUGH, oaMessages_)<<cut2rowIdx[cutIdx]<<violation<<CoinMessageEol;
        continue;}
        if(oaHandler_->logLevel() > 0)
          oaHandler_->message(VIOLATED_OA_CUT_GENERATED, oaMessages_)<<cut2rowIdx[cutIdx]<<violation<<CoinMessageEol;
    }
    else if (oaHandler_->logLevel() > 0)
      oaHandler_->message(OA_CUT_GENERATED, oaMessages_)<<cut2rowIdx[cutIdx]<<CoinMessageEol;
  OsiRowCut newCut;
    //    if(lb[i]>-1e20) assert (ub[i]>1e20);

    if (IsValid(cutStrengthener_)) {
      const int& rowIdx = cut2rowIdx[cutIdx];
      bool retval =
	cutStrengthener_->ComputeCuts(cs, GetRawPtr(tminlp_),
				       GetRawPtr(problem_), rowIdx,
				       cuts[cutIdx], lb[cutIdx], ub[cutIdx], g[rowIdx],
				       rowLower[rowIdx], rowUpper[rowIdx],
				       n, x, infty);
      if (!retval) {
	(*messageHandler()) << "error in cutStrengthener_->ComputeCuts\n";
	//exit(-2);
      }
    }
    if(global) {
      newCut.setGloballyValidAsInteger(1);
    }
    newCut.setEffectiveness(99.99e99);
    newCut.setLb(lb[cutIdx]);
    newCut.setUb(ub[cutIdx]);
    newCut.setRow(cuts[cutIdx]);
    //    CftValidator validator;
    //    validator(newCut);
    if(oaHandler_->logLevel()>2){
      oaHandler_->print(newCut);}
    cs.insert(newCut);
  }

  delete[] g;
  delete [] cuts;
  delete [] row2cutIdx;
  delete [] cut2rowIdx;

  if(getObj) { // Get the objective cuts
    double * obj = new double [n];
    tminlp_->eval_grad_f(n, x, 1,obj);
    double f;
    tminlp_->eval_f(n, x, 1, f);

    CoinPackedVector v;
    v.reserve(n);
    lb[nNonLinear_] = -f;
    ub[nNonLinear_] = -f;
    //double minCoeff = 1e50;
    for(int i = 0; i<n ; i++) {
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
    //Compute cut violation
    bool genCut = true;
    if(x2 != NULL) {
      double rhs = v.dotProduct(x2);
      double violation = max(0., rhs - ub[nNonLinear_]);
      if(violation < theta) genCut = false;
    }
    if(genCut) {
      if (IsValid(cutStrengthener_)) {
	lb[nNonLinear_] = -infty;
	bool retval =
	  cutStrengthener_->ComputeCuts(cs, GetRawPtr(tminlp_),
					 GetRawPtr(problem_), -1,
					 v, lb[nNonLinear_], ub[nNonLinear_],
					 ub[nNonLinear_], -infty, 0.,
					 n, x, infty);
	if (!retval) {
	  std::cerr << "error in cutStrengthener_->ComputeCuts\n";
	  //exit(-2);
	}
      }
      OsiRowCut newCut;
      if(global)
	newCut.setGloballyValidAsInteger(1);
      newCut.setEffectiveness(99.99e99);
      newCut.setRow(v);
      newCut.setLb(-COIN_DBL_MAX/*Infinity*/);
      newCut.setUb(ub[nNonLinear_]);
      //     CftValidator validator;
      //     validator(newCut);
      cs.insert(newCut);
    }
    delete [] obj;
    }

  delete []lb;
  delete[]ub;
}



/** Get the outer approximation of a single constraint at the point x.
*/
void
OsiTMINLPInterface::getConstraintOuterApproximation(OsiCuts &cs, int rowIdx, 
                                                    const double * x, 
                                                    const double * x2, bool global)
{
  double g;
  int * indices = new int[getNumCols()];
  double * values = new double[getNumCols()];
  int nnz;
  tminlp_->eval_grad_gi(getNumCols(), x, 1, rowIdx, nnz, indices, values);
  tminlp_->eval_gi(getNumCols(),x,1, rowIdx, g);

  CoinPackedVector cut;
  double lb;
  double ub;


  const double rowLower = getRowLower()[rowIdx];
  const double rowUpper = getRowUpper()[rowIdx];
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double dual = (getRowPrice() + 2 * getNumCols())[rowIdx];
  double infty = getInfinity();
  double nlp_infty = infty_;
  
  if(rowLower > - nlp_infty)
    lb = rowLower - g;
  else
    lb = - infty;
  if(rowUpper < nlp_infty)
    ub = rowUpper - g;
  else
    ub = infty;
  if(rowLower > -infty && rowUpper < infty)
  {
    if(dual >= 0)// <= inequality
      lb = - infty;
    if(dual <= 0)// >= inequality
      ub = infty;
  }

  for(int i = 0 ; i < nnz; i++) {
     const int &colIdx = indices[i];
      //"clean" coefficient
      if(cleanNnz(values[i],colLower[colIdx], colUpper[colIdx],
		  rowLower, rowUpper,
		  x[colIdx],
		  lb,
		  ub, tiny_, veryTiny_)) {
        cut.insert(colIdx,values[i]);
        if(lb > - infty)
          lb += values[i] * x[colIdx];
        if(ub < infty)
	  ub += values[i] * x[colIdx];
    }
  }

  OsiRowCut newCut;

  if(global) {
    newCut.setGloballyValidAsInteger(1);
  }
  newCut.setEffectiveness(99.99e99);
  newCut.setLb(lb);
  newCut.setUb(ub);
  newCut.setRow(cut);
  cs.insert(newCut);

  delete [] indices;
  delete [] values;
}


double
OsiTMINLPInterface::getFeasibilityOuterApproximation(int n,const double * x_bar,const int *inds, OsiCuts &cs, bool addOnlyViolated, bool global)
{
  if(! IsValid(feasibilityProblem_)) {
    throw SimpleError("No feasibility problem","getFeasibilityOuterApproximation");
  }
  feasibilityProblem_->set_dist2point_obj(n,(const Number *) x_bar,(const Index *) inds);
  nCallOptimizeTNLP_++;
  totalNlpSolveTime_-=CoinCpuTime();
  SmartPtr<TNLPSolver> app2 = app_->clone();
  app2->options()->SetIntegerValue("print_level", (Index) 0);
  optimizationStatus_ = app2->OptimizeTNLP(GetRawPtr(feasibilityProblem_));
  totalNlpSolveTime_+=CoinCpuTime();
  getOuterApproximation(cs, getColSolution(), 0, (addOnlyViolated? x_bar:NULL)
			, global);
  hasBeenOptimized_=true;
  return getObjValue();
}


static bool WarnedForNonConvexOa=false;


void
OsiTMINLPInterface::extractLinearRelaxation(OsiSolverInterface &si, bool getObj,
                                            bool solveNlp)
{
  if(solveNlp)
    initialSolve();

  double * rowLow = NULL;
  double * rowUp = NULL;

  int n;
  int m;
  int nnz_jac_g;
  int nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
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
  int * nonBindings = new int[m];//store non binding constraints (which are to be removed from OA)
  int numNonBindings = 0;
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  const double * colLower = getColLower();
  const double * colUpper = getColUpper();
  const double * duals = getRowPrice() + 2*n;
  assert(m==getNumRows());
  double infty = getInfinity();
  double nlp_infty = infty_;
  
  for(int i = 0 ; i < m ; i++) {
    if(constTypes_[i] == TNLP::NON_LINEAR) {
      //If constraint is range not binding prepare to remove it
      if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty && fabs(duals[i]) == 0.)
      {
        nonBindings[numNonBindings++] = i;
        continue;
      }
      else
        if(rowLower[i] > - nlp_infty){
          rowLow[i] = (rowLower[i] - g[i]) - 1e-07;
          if(! WarnedForNonConvexOa && rowUpper[i] < nlp_infty){
             messageHandler()->message(WARNING_NON_CONVEX_OA, messages_)<<CoinMessageEol;
             WarnedForNonConvexOa = true;
          }
        }
      else
        rowLow[i] = - infty;
      if(rowUpper[i] < nlp_infty)
        rowUp[i] =  (rowUpper[i] - g[i]) + 1e-07;
      else
        rowUp[i] = infty;
      
      //If equality or ranged constraint only add one side by looking at sign of dual multiplier
      if(rowLower[i] > -nlp_infty && rowUpper[i] < nlp_infty)
      {
        if(duals[i] >= 0.)// <= inequality
          rowLow[i] = - infty;
        if(duals[i] <= 0.)// >= inequality
          rowUp[i] = infty;
      }
    }
    else {
      rowLow[i] = (rowLower[i] - g[i]);
      rowUp[i] =  (rowUpper[i] - g[i]);
    }
  }

  
  
  //Then convert everything to a CoinPackedMatrix
  //Go through values, clean coefficients and fix bounds
  for(int i = 0 ; i < nnz_jac_g ; i++) {
    if(constTypes_[jRow_[i]] == TNLP::LINEAR //Always accept coefficients from linear constraints
       || //For other clean tinys
       cleanNnz(jValues_[i],colLower[jCol_[i]], colUpper[jCol_[i]],
                rowLower[jRow_[i]], rowUpper[jRow_[i]],
                getColSolution()[jCol_[i]],
                rowLow[jRow_[i]],
                rowUp[jRow_[i]], tiny_, veryTiny_)) {      
      rowLow[jRow_[i]] += jValues_[i] * getColSolution()[jCol_ [i]];
      rowUp[jRow_[i]] += jValues_[i] *getColSolution()[jCol_[i]];
    }
  }
  CoinPackedMatrix mat(true, jRow_, jCol_, jValues_, nnz_jac_g);
  mat.setDimensions(m,n); // In case matrix was empty, this should be enough
  

  //remove non-bindings equality constraints
  mat.deleteRows(numNonBindings, nonBindings);
  
  int numcols=getNumCols();
  double *obj = new double[numcols];
  for(int i = 0 ; i < numcols ; i++)
    obj[i] = 0.;
  
  
  si.loadProblem(mat, getColLower(), getColUpper(), obj, rowLow, rowUp);
  delete [] rowLow;
  delete [] rowUp;
  delete [] nonBindings;
  delete [] g;
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
    optimizationStatus_ = app_->ReOptimizeTNLP(GetRawPtr(problem_));
  else
    optimizationStatus_ = app_->OptimizeTNLP(GetRawPtr(problem_));
  totalNlpSolveTime_+=CoinCpuTime();
  nCallOptimizeTNLP_++;
  hasBeenOptimized_ = true;
  
  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    //app_->Options()->SetIntegerValue("print_level",0, true, true);
    app_->options()->SetStringValue("print_user_options","no", false, true);
  }
  
  
#if 1
  if(optimizationStatus_ == TNLPSolver::notEnoughFreedom)//Too few degrees of freedom
  {
    (*messageHandler())<<"Too few degrees of freedom...."<<CoinMessageEol;
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
    if(numcols - numberFixed > numberEqualities || numcols < numberEqualities)
    {
      std::string probName;
      getStrParam(OsiProbName, probName);
      throw newUnsolvedError(app_->errorCode(), problem_, probName);
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
    if(!app_->isRecoverable(optimizationStatus_))//Solver failed and the error can not be recovered, throw it
    {
      std::string probName;
      getStrParam(OsiProbName, probName);
      throw newUnsolvedError(app_->errorCode(), problem_, probName);
    }
  try{
    totalIterations_ += app_->IterationCount();
  }
  catch(SimpleError &E)
  {
    if (throwOnFailure)//something failed throw
    {
      throw SimpleError("No statistics available from Ipopt",whereFrom);
    }
    else {
      return;
    }
  }
  if(problem_->hasUpperBoundingObjective()){//Check if solution is integer and recompute objective value using alternative objective function
    const double * sol = getColSolution();
    bool integerSol = true;
    double intTol = 1e-08;
    if(objects()){
      int nObjects = numberObjects();
      OsiObject ** object = objects();
      for(int i = 0 ; i< nObjects ; i++){
        int dummy;
        if(object[i]->infeasibility(this,dummy) > intTol)
        {
          integerSol=false;
          break;
        }
      }
    }
    else{//Only works for integer constraints
      int numcols = getNumCols();
      for(int i = 0 ; i < numcols ; i++){
        if(isInteger(i) || isBinary(i)){
          if(fabs(sol[i] - floor(sol[i]+0.5)) > intTol){
            integerSol = false;
            break;
          }
        }
      }
    }
    if(integerSol){
      problem_->evaluateUpperBoundingFunction(sol);
      messageHandler()->message(ALTERNATE_OBJECTIVE, messages_)
      <<getObjValue()<<CoinMessageEol;
    }
  }
  messageHandler()->message(IPOPT_SUMMARY, messages_)
    <<whereFrom<<optimizationStatus_<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;
  
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
    app_->options()->GetEnumValue("print_user_options",printOptions,"bonmin.");
    if(printOptions)
      app_->options()->SetStringValue("print_user_options","yes");
  }
  app_->disableWarmStart(); 
  solveAndCheckErrors(0,1,"initialSolve");
  
  //Options should have been printed if not done already turn off Ipopt output
  if(!hasPrintedOptions) {
    hasPrintedOptions = 1;
    app_->options()->SetStringValue("print_user_options","no");
    app_->options()->SetIntegerValue("print_level",0);
  }
  
  const char * status=OPT_SYMB;
  ;
  if(isAbandoned()) status=FAILED_SYMB;
  else if(isProvenPrimalInfeasible()) status=INFEAS_SYMB;
  messageHandler()->message(LOG_FIRST_LINE, messages_)<<nCallOptimizeTNLP_
						      <<status<<getObjValue()<<app_->IterationCount()<<app_->CPUTime()<<CoinMessageEol;
  
  int numRetry = firstSolve_ ? numRetryInitial_ : numRetryResolve_;
  if(isAbandoned()) {
    resolveForRobustness(numRetryUnsolved_);
  }
  else if(numRetry)
    {
      resolveForCost(numRetry);
      /** Only do it once at the root.*/
      numRetryInitial_ = 0;
    }
  firstSolve_ = false;
}

/** Resolve the continuous relaxation after problem modification.
 * \note for Ipopt, same as resolve */
void
OsiTMINLPInterface::resolve()
{
  assert(IsValid(app_));
  assert(IsValid(problem_));
  if (INT_BIAS > 0.) {
    app_->options()->SetStringValue("warm_start_same_structure", "yes");
  }
  else {
    app_->options()->SetStringValue("warm_start_same_structure", "no");
  }

  if(problem_->duals_init() != NULL)
    app_->enableWarmStart();
  else app_->disableWarmStart();
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
  else if(numRetryResolve_ ||
	  (numRetryInfeasibles_ && isProvenPrimalInfeasible() ))
	  resolveForCost(max(numRetryResolve_, numRetryInfeasibles_));
  
}


////////////////////////////////////////////////////////////////
// Methods returning info on how the solution process terminated  //
////////////////////////////////////////////////////////////////
/// Are there a numerical difficulties?
bool OsiTMINLPInterface::isAbandoned() const
{
  return (
        (optimizationStatus_==TNLPSolver::iterationLimit)||
        (optimizationStatus_==TNLPSolver::computationError)||
        (optimizationStatus_==TNLPSolver::illDefinedProblem)||
        (optimizationStatus_==TNLPSolver::illegalOption)||
        (optimizationStatus_==TNLPSolver::externalException)||
        (optimizationStatus_==TNLPSolver::exception)
      );
}

/// Is optimality proven?
bool OsiTMINLPInterface::isProvenOptimal() const
{
  return (optimizationStatus_==TNLPSolver::solvedOptimal ||
	  optimizationStatus_==TNLPSolver::solvedOptimalTol);
}
/// Is primal infeasiblity proven?
bool OsiTMINLPInterface::isProvenPrimalInfeasible() const
{
  return (optimizationStatus_ == TNLPSolver::provenInfeasible);
}
/// Is dual infeasiblity proven?
bool OsiTMINLPInterface::isProvenDualInfeasible() const
{
  return (optimizationStatus_ == TNLPSolver::unbounded);
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
  return (optimizationStatus_==TNLPSolver::unbounded);

}
/// Iteration limit reached?
bool OsiTMINLPInterface::isIterationLimitReached() const
{
  return (optimizationStatus_==TNLPSolver::iterationLimit);
}

void
OsiTMINLPInterface::extractInterfaceParams()
{
  if (IsValid(app_)) {
    int logLevel;
    app_->options()->GetIntegerValue("nlp_log_level", logLevel,"bonmin.");
    messageHandler()->setLogLevel(logLevel);

#ifdef COIN_HAS_FILTERSQP
    FilterSolver * filter = dynamic_cast<FilterSolver *>(GetRawPtr(app_));

    bool is_given =
#endif
      app_->options()->GetNumericValue("max_random_point_radius",maxRandomRadius_,"bonmin.");

#ifdef COIN_HAS_FILTERSQP
    if(filter && !is_given){
      // overwriting default for filter
      maxRandomRadius_ = 10.;
    }
#endif
   
   int oaCgLogLevel = 0;
   app_->options()->GetIntegerValue("oa_cuts_log_level", oaCgLogLevel,"bonmin.");
   oaHandler_->setLogLevel(oaCgLogLevel); 
    
    int exposeWs = false;
    app_->options()->GetEnumValue("warm_start", exposeWs, "bonmin.");
    setExposeWarmStart(exposeWs > 0);
    
    app_->options()->GetIntegerValue("num_retry_unsolved_random_point", numRetryUnsolved_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_root", numRetryInitial_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_node", numRetryResolve_,"bonmin.");
    app_->options()->GetIntegerValue("num_resolve_at_infeasibles", numRetryInfeasibles_,"bonmin.");
    app_->options()->GetIntegerValue("num_iterations_suspect", numIterationSuspect_,"bonmin.");
    app_->options()->GetEnumValue("nlp_failure_behavior",pretendFailIsInfeasible_,"bonmin.");
    app_->options()->GetNumericValue
    ("warm_start_bound_frac" ,pushValue_,"bonmin.");
    app_->options()->GetNumericValue("tiny_element",tiny_,"bonmin.");
    app_->options()->GetNumericValue("very_tiny_element",veryTiny_,"bonmin.");
    app_->options()->GetNumericValue("random_point_perturbation_interval",max_perturbation_,"bonmin.");
    app_->options()->GetEnumValue("random_point_type",randomGenerationType_,"bonmin.");
    int cut_strengthening_type;
    app_->options()->GetEnumValue("cut_strengthening_type", cut_strengthening_type,"bonmin.");

    if (cut_strengthening_type != CS_None) {
      // TNLP solver to be used in the cut strengthener
      cutStrengthener_ = new CutStrengthener(app_->clone(), app_->options());
    }
  }
}

void
OsiTMINLPInterface::SetStrongBrachingSolver(SmartPtr<StrongBranchingSolver> strong_branching_solver)
{
  strong_branching_solver_ = strong_branching_solver;
}

  //#define STRONG_COMPARE
#ifdef STRONG_COMPARE
  static double objorig;
#endif

void
OsiTMINLPInterface::markHotStart()
{
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::markHotStart();
    objorig = getObjValue();
#endif
    optimizationStatusBeforeHotStart_ = optimizationStatus_;
    strong_branching_solver_->markHotStart(this);
  }
  else {
    // Default Implementation
    OsiSolverInterface::markHotStart();
  }
}

void
OsiTMINLPInterface::solveFromHotStart()
{
#if 0
  printf("========= 1111111111111 ==============\n");
  for (int i=0; i<getNumCols(); i++) {
    printf("xL[%3d] = %15.8e  xU[%3d] = %15.8e\n", i, getColLower()[i], i, getColUpper()[i]);
  }
#endif
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::solveFromHotStart();
    double obj_nlp = getObjValue() - objorig;
#endif
    optimizationStatus_ = strong_branching_solver_->solveFromHotStart(this);
#ifdef STRONG_COMPARE
    double obj_other = getObjValue() - objorig;
    printf("AWDEBUG: Strong Branching results: NLP = %15.8e Other = %15.8e\n",
	   obj_nlp, obj_other);
#endif
  }
  else {
    // Default Implementation
    OsiSolverInterface::solveFromHotStart();
  }
}

void
OsiTMINLPInterface::unmarkHotStart()
{
  if (IsValid(strong_branching_solver_)) {
#ifdef STRONG_COMPARE
    // AWDEBUG
    OsiSolverInterface::unmarkHotStart();
#endif
    strong_branching_solver_->unmarkHotStart(this);
    optimizationStatus_ = optimizationStatusBeforeHotStart_;
  }
  else {
    // Default Implementation
    OsiSolverInterface::unmarkHotStart();
  }
}

const double * OsiTMINLPInterface::getObjCoefficients() const
{
  const int n = getNumCols();
  delete [] obj_;
  obj_ = NULL;
  obj_ = new double[n];

  bool new_x = true;
  const double* x_sol = problem_->x_sol();
  bool retval = problem_->eval_grad_f(n, x_sol, new_x, obj_);
  
  if (!retval) {
    // Let's see if that happens - it will cause a crash
    printf("ERROR WHILE EVALUATING GRAD_F in OsiTMINLPInterface::getObjCoefficients()\n");
    delete [] obj_;
    obj_ = NULL;
  }

  return obj_;
}


}/** end namespace Bonmin*/

