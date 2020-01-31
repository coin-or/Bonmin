// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007

#include "BonminConfig.h"
#ifdef COIN_HAS_FILTERSQP
# include "BonFilterSolver.hpp"
#endif
#include "BonBabSetupBase.hpp"
#include <climits>
#include <fstream>
#include <sstream>

#include "BonDiver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"
#include "BonChooseVariable.hpp"
#include "BonTMINLP2Quad.hpp"
#include "BonTMINLPLinObj.hpp"
namespace Bonmin
{
  int BabSetupBase::defaultIntParam_[BabSetupBase::NumberIntParam] = {
        1 /* BabLogLevel*/,
        100 /* BabLogInterval*/,
        2 /* MaxFailures.*/,
        0 /* FailureBehavior.*/,
        0 /* MaxInfeasible*/,
        5 /*NumberStrong*/,
        2 /* MinReliability*/,
        COIN_INT_MAX /* MaxNodes*/,
        COIN_INT_MAX /* MaxSolutions*/,
        COIN_INT_MAX /* MaxIterations*/,
        0 /* SpecialOption*/,
        0 /* DisableSos.*/,
        1 /* numCutPasses.*/,
        20 /* numCutPassesAtRoot.*/,
        0 /* log level at root.*/
      };


  double BabSetupBase::defaultDoubleParam_[BabSetupBase::NumberDoubleParam] = {
        0 /* CutoffDecr*/,
        COIN_DBL_MAX /* Cutoff*/,
        0 /* AllowableGap*/,
        0 /*AllowableFractionGap*/,
        1e-09 /*IntTol*/,
        COIN_DBL_MAX /* MaxTime*/,
      };

  BabSetupBase::BabSetupBase(const CoinMessageHandler * handler):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(NULL),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(),
      treeTraversalMethod_(),
      objects_(0),
      journalist_(NULL),
      options_(NULL),
      roptions_(NULL),
      readOptions_(false),
      messageHandler_(NULL),
      prefix_("bonmin.")
  {
    CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
    CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
    if(handler) messageHandler_ = handler->clone();
  }

  /** Copy constructor. */
  BabSetupBase::BabSetupBase(const BabSetupBase & other):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(other.linearizer_),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(other.nodeComparisonMethod_),
      treeTraversalMethod_(other.treeTraversalMethod_),
      objects_(other.objects_),
      journalist_(other.journalist_),
      options_(NULL),
      roptions_(other.roptions_),
      readOptions_(other.readOptions_),
      messageHandler_(NULL),
      prefix_(other.prefix_)
  {
    if (other.nonlinearSolver_) {
      nonlinearSolver_ = static_cast<OsiTMINLPInterface *>(other.nonlinearSolver_->clone());
    }
    if (other.continuousSolver_) {
      continuousSolver_ = other.continuousSolver_->clone();
    }
    if (other.messageHandler_) {
      messageHandler_ = other.messageHandler_->clone();
      continuousSolver_->passInMessageHandler(messageHandler_);
    }
    for (CuttingMethods::const_iterator i = other.cutGenerators_.begin() ; i != other.cutGenerators_.end() ; i++) {
      cutGenerators_.push_back(*i);
      cutGenerators_.back().cgl = cutGenerators_.back().cgl->clone();
    }

    for (HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++) {
      heuristics_.push_back(*i);
      heuristics_.back().heuristic = i->heuristic->clone();
    }
  
    if(other.branchingMethod_ != NULL)
      branchingMethod_ = other.branchingMethod_->clone();
    if (IsValid(other.options_)) {
      options_ = new Ipopt::OptionsList;
      *options_ = *other.options_;
    }
    CoinCopyN(other.intParam_, NumberIntParam, intParam_);
    CoinCopyN(other.doubleParam_, NumberDoubleParam, doubleParam_);
    for (unsigned int i = 0 ; i < objects_.size() ; i++) {
      objects_[i]->clone();
    }
  }

  /** Copy constructor with change of nlp. */
  BabSetupBase::BabSetupBase(const BabSetupBase & other,
                             OsiTMINLPInterface &nlp):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(other.linearizer_),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(other.nodeComparisonMethod_),
      treeTraversalMethod_(other.treeTraversalMethod_),
      objects_(other.objects_),
      journalist_(other.journalist_),
      options_(NULL),
      roptions_(other.roptions_),
      readOptions_(other.readOptions_),
      messageHandler_(NULL),
      prefix_(other.prefix_)
  {
    nonlinearSolver_ = &nlp;
    if (other.continuousSolver_ != other.nonlinearSolver_) {
      continuousSolver_ = NULL;
    }
    else
      continuousSolver_ = nonlinearSolver_;
    if (other.messageHandler_) {
      messageHandler_ = other.messageHandler_->clone();
      continuousSolver_->passInMessageHandler(messageHandler_);
    }
    for (CuttingMethods::const_iterator i = other.cutGenerators_.begin() ; i != other.cutGenerators_.end() ; i++) {
      cutGenerators_.push_back(*i);
      cutGenerators_.back().cgl = cutGenerators_.back().cgl->clone();
    }

    for (HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++) {
      heuristics_.push_back(*i);
      heuristics_.back().heuristic = i->heuristic->clone();
    }
  
    if(other.branchingMethod_ != NULL)
      branchingMethod_ = other.branchingMethod_->clone();
    if (IsValid(other.options_)) {
      options_ = new Ipopt::OptionsList;
      *options_ = *other.options_;
    }
    CoinCopyN(other.intParam_, NumberIntParam, intParam_);
    CoinCopyN(other.doubleParam_, NumberDoubleParam, doubleParam_);
    for (unsigned int i = 0 ; i < objects_.size() ; i++) {
      objects_[i]->clone();
    }
  }
  /** Copy constructor with change of nlp. */
  BabSetupBase::BabSetupBase(const BabSetupBase & other,
                             OsiTMINLPInterface &nlp,
                             const std::string & prefix):
      nonlinearSolver_(other.nonlinearSolver_),
      continuousSolver_(NULL),
      linearizer_(other.linearizer_),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(),
      treeTraversalMethod_(),
      objects_(other.objects_),
      journalist_(other.journalist_),
      options_(NULL),
      roptions_(other.roptions_),
      readOptions_(other.readOptions_),
      messageHandler_(NULL),
      prefix_(prefix)
  {
    nonlinearSolver_ = &nlp;
    if (IsValid(other.options_)) {
      options_ = new Ipopt::OptionsList;
      *options_ = *other.options_;
    }
    if (other.messageHandler_) {
      messageHandler_ = other.messageHandler_->clone();
      nonlinearSolver_->passInMessageHandler(messageHandler_);
    }
    CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
    CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
    gatherParametersValues(options_);
    for (unsigned int i = 0 ; i < objects_.size() ; i++) {
      objects_[i]->clone();
    }
  }

  BabSetupBase::BabSetupBase(Ipopt::SmartPtr<TMINLP> tminlp, const CoinMessageHandler * handler):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(NULL),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(),
      treeTraversalMethod_(),
      objects_(0),
      readOptions_(false),
      messageHandler_(NULL),
      prefix_("bonmin.")
  {
    CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
    CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
    if(handler) messageHandler_ = handler->clone();
    use(tminlp);
  }


  /** Make a copy with solver replace by one passed .*/
  BabSetupBase *
  BabSetupBase::clone(OsiTMINLPInterface&nlp)const {
     throw(CoinError("BabSetupBase", "CloneWithOtherNlp","Not implemented"));
  }


  void
  BabSetupBase::use(Ipopt::SmartPtr<TMINLP> tminlp)
  {
    readOptionsFile();
    assert(IsValid(tminlp));
    nonlinearSolver_ = new OsiTMINLPInterface;
    int ival;
    options_->GetEnumValue("enable_dynamic_nlp", ival, "bonmin.");
    if(ival && ! tminlp->hasLinearObjective()){
      Ipopt::SmartPtr<Bonmin::TMINLPLinObj> linObj =
                        new Bonmin::TMINLPLinObj;
      linObj->setTminlp(GetRawPtr(tminlp));
      tminlp = GetRawPtr(linObj);
    }
    nonlinearSolver_->initialize(roptions_, options_, journalist_, prefix(), tminlp);
    if(messageHandler_ != NULL)
      nonlinearSolver_->passInMessageHandler(messageHandler_);
    else
      messageHandler_ = nonlinearSolver_->messageHandler()->clone();
    if (ival){
      nonlinearSolver_->use(new Bonmin::TMINLP2TNLPQuadCuts(tminlp));
    }
  }

  BabSetupBase::BabSetupBase(const OsiTMINLPInterface& nlp):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(NULL),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(),
      treeTraversalMethod_(),
      objects_(0),
      journalist_(NULL),
      options_(NULL),
      roptions_(NULL),
      readOptions_(false),
      messageHandler_(NULL),
      prefix_("bonmin.")
  {
    CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
    CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
    use(nlp);
  }

  void
  BabSetupBase::use(const OsiTMINLPInterface& nlp)
  {
    nonlinearSolver_ = dynamic_cast<OsiTMINLPInterface *>(nlp.clone());
    options_ = nonlinearSolver_->solver()->options();
    roptions_ = nonlinearSolver_->solver()->roptions();
    journalist_ = nonlinearSolver_->solver()->journalist();
    if(messageHandler_ != NULL ) delete messageHandler_;		
    messageHandler_ = nlp.messageHandler()->clone();
    readOptions_ = true;
  }

  BabSetupBase::BabSetupBase( Ipopt::SmartPtr<TNLPSolver> app):
      nonlinearSolver_(NULL),
      continuousSolver_(NULL),
      linearizer_(NULL),
      cutGenerators_(),
      heuristics_(),
      branchingMethod_(NULL),
      nodeComparisonMethod_(),
      treeTraversalMethod_(),
      objects_(0),
      journalist_(app->journalist()),
      options_(app->options()),
      roptions_(app->roptions()),
      readOptions_(true),
      messageHandler_(NULL),
      prefix_("bonmin.")
  {
    CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
    CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
  }


  BabSetupBase::~BabSetupBase()
  {
    if (nonlinearSolver_ != continuousSolver_) {
      delete nonlinearSolver_;
    }
    delete continuousSolver_;
    delete branchingMethod_;
    for (CuttingMethods::iterator i = cutGenerators_.begin() ; i != cutGenerators_.end() ; i++) {
      delete i->cgl;
      i->cgl = NULL;
    }

    for (HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++) {
      delete i->heuristic;
    }

    for (unsigned int i = 0 ; i < objects_.size() ; i++) {
      delete objects_[i];
    }

    if(messageHandler_)
      delete messageHandler_;
  }


  void
  BabSetupBase::gatherParametersValues(Ipopt::SmartPtr<Ipopt::OptionsList> options)
  {

    options->GetIntegerValue("bb_log_level",intParam_[BabLogLevel],prefix_.c_str());
    options->GetIntegerValue("bb_log_interval",intParam_[BabLogInterval],prefix_.c_str());
    options->GetIntegerValue("max_consecutive_failures",intParam_[MaxFailures],prefix_.c_str());
    options->GetEnumValue("nlp_failure_behavior",intParam_[FailureBehavior],prefix_.c_str());
    options->GetIntegerValue("max_consecutive_infeasible",intParam_[MaxInfeasible],prefix_.c_str());
    options->GetIntegerValue("number_strong_branch",intParam_[NumberStrong],prefix_.c_str());
    options->GetIntegerValue("number_before_trust",intParam_[MinReliability],prefix_.c_str());
    options->GetIntegerValue("node_limit",intParam_[MaxNodes],prefix_.c_str());
    options->GetIntegerValue("solution_limit",intParam_[MaxSolutions],prefix_.c_str());
    options->GetIntegerValue("iteration_limit",intParam_[MaxIterations],prefix_.c_str());
    options->GetEnumValue("sos_constraints",intParam_[DisableSos],prefix_.c_str());
    options->GetIntegerValue("num_cut_passes",intParam_[NumCutPasses],prefix_.c_str());
    options->GetIntegerValue("num_cut_passes_at_root",intParam_[NumCutPassesAtRoot],prefix_.c_str());
    options->GetIntegerValue("nlp_log_at_root",intParam_[RootLogLevel],prefix_.c_str());

    options->GetNumericValue("cutoff_decr",doubleParam_[CutoffDecr],prefix_.c_str());
    options->GetNumericValue("cutoff",doubleParam_[Cutoff],prefix_.c_str());
    options->GetNumericValue("allowable_gap",doubleParam_[AllowableGap],prefix_.c_str());
    options->GetNumericValue("allowable_fraction_gap",doubleParam_[AllowableFractionGap],prefix_.c_str());
    options->GetNumericValue("integer_tolerance",doubleParam_[IntTol],prefix_.c_str());
    options->GetNumericValue("time_limit", doubleParam_[MaxTime],prefix_.c_str());

    int ival;
    int seed = 0;
    ival = options->GetIntegerValue("random_generator_seed",seed,prefix_.c_str());
    if(seed == -1)
      CoinSeedRandom((int)CoinGetTimeOfDay());
    else if (ival != 0) CoinSeedRandom(seed);

    options->GetEnumValue("node_comparison",ival,prefix_.c_str());
    nodeComparisonMethod_ = NodeComparison(ival);

    options->GetEnumValue("tree_search_strategy", ival, prefix_.c_str());
    treeTraversalMethod_ = TreeTraversal(ival);

    int varSelection;
    options->GetEnumValue("variable_selection",varSelection,prefix_.c_str());
    // Set branching strategy
    if (varSelection == MOST_FRACTIONAL) {
      intParam_[NumberStrong] = 0;
      intParam_[MinReliability] = 0;
      options_->SetIntegerValue("bonmin.number_strong_branch",intParam_[BabSetupBase::NumberStrong],true,true);
      options_->SetIntegerValue("bonmin.number_before_trust",intParam_[BabSetupBase::MinReliability],true,true);
    }
    else if (varSelection == RELIABILITY_BRANCHING) {
      intParam_[MinReliability] = 10;
      options_->SetIntegerValue("bonmin.number_before_trust",intParam_[BabSetupBase::MinReliability],true,true);
    }
  }

  void
  BabSetupBase::registerOptions()
  {
    registerAllOptions(roptions_);
  }

  void
  BabSetupBase::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    OsiTMINLPInterface::registerOptions(roptions);
    /* BabSetup options.*/
    roptions->SetRegisteringCategory("Output and Loglevel", RegisteredOptions::BonminCategory);

    roptions->AddBoundedIntegerOption("bb_log_level",
        "specify main branch-and-bound log level.",
        0,5,1,
        "Set the level of output of the branch-and-bound : "
        "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                     );
    roptions->setOptionExtraInfo("bb_log_level", 127);

    roptions->AddLowerBoundedIntegerOption("bb_log_interval",
        "Interval at which node level output is printed.",
        0,100,
        "Set the interval (in terms of number of nodes) at which "
        "a log on node resolutions (consisting of lower and upper bounds) is given.");
    roptions->setOptionExtraInfo("bb_log_interval", 127);

    roptions->AddBoundedIntegerOption("lp_log_level",
        "specify LP log level.",
        0,4,0,
        "Set the level of output of the linear programming sub-solver in B-Hyb or B-QG : "
        "0 - none, 1 - minimal, 2 - normal low, 3 - normal high, 4 - verbose"
                                     );
    roptions->setOptionExtraInfo("lp_log_level", 119);

    roptions->AddBoundedIntegerOption("nlp_log_at_root","specify a different log level "
                                           "for root relaxation.",
                                            0,12,0,
                                            "");
    roptions->setOptionExtraInfo("nlp_log_at_root",63);

    roptions->SetRegisteringCategory("Branch-and-bound options", RegisteredOptions::BonminCategory);

  roptions->AddLowerBoundedIntegerOption
  ("random_generator_seed",
   "Set seed for random number generator (a value of -1 sets seeds to time since Epoch).",
   -1,0,
   "");
  roptions->setOptionExtraInfo("random_generator_seed",127);

    roptions->AddLowerBoundedNumberOption("time_limit",
        "Set the global maximum computation time (in secs) for the algorithm.",
        0.,0,1e10,
        "");
    roptions->setOptionExtraInfo("time_limit", 127);

    roptions->AddLowerBoundedIntegerOption("node_limit",
        "Set the maximum number of nodes explored in the branch-and-bound search.",
        0,COIN_INT_MAX,
        "");
    roptions->setOptionExtraInfo("node_limit", 127);

    roptions->AddLowerBoundedIntegerOption("iteration_limit",
        "Set the cumulative maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound.",
        0,COIN_INT_MAX,
        "value 0 deactivates option.");
    roptions->setOptionExtraInfo("iteration_limit", 127);

    roptions->AddLowerBoundedIntegerOption("solution_limit",
        "Abort after that much integer feasible solution have been found by algorithm",
        0,COIN_INT_MAX,
        "value 0 deactivates option");
    roptions->setOptionExtraInfo("solution_limit", 127);

    roptions->AddLowerBoundedNumberOption("integer_tolerance",
        "Set integer tolerance.",
        0.,1,1e-06,
        "Any number within that value of an integer is considered integer.");
    roptions->setOptionExtraInfo("integer_tolerance", 127);

    roptions->AddBoundedNumberOption("allowable_gap",
        "Specify the value of absolute gap under which the algorithm stops.",
        -1.e20,0,1.e20,0,0.,
        "Stop the tree search when the gap between the objective value of the best known solution"
        " and the best bound on the objective of any solution is less than this.");
    roptions->setOptionExtraInfo("allowable_gap", 127);

    roptions->AddBoundedNumberOption("allowable_fraction_gap",
        "Specify the value of relative gap under which the algorithm stops.",
        -1.e20,0,1.e20,0,0.0,
        "Stop the tree search when the gap between the objective value of the best known solution "
        "and the best bound on the objective of any solution is less than this "
        "fraction of the absolute value of the best known solution value.");
    roptions->setOptionExtraInfo("allowable_fraction_gap", 127);

    roptions->AddBoundedNumberOption("cutoff",
        "Specify cutoff value.",
        -1e100,0,1e100,0,1e100,
        "cutoff should be the value of a feasible solution known by the user "
        "(if any). The algorithm will only look for solutions better than cutoff.");
    roptions->setOptionExtraInfo("cutoff", 127);


    roptions->AddBoundedNumberOption("cutoff_decr",
        "Specify cutoff decrement.",
        -1e10,0,1e10,0,1e-05,
        "Specify the amount by which cutoff is decremented below "
        "a new best upper-bound"
        " (usually a small positive value but in non-convex problems it may be a negative value).");
    roptions->setOptionExtraInfo("cutoff_decr", 127);


    roptions->AddStringOption5("node_comparison",
        "Choose the node selection strategy.",
        "best-bound",
        "best-bound", "choose node with the smallest bound,",
        "depth-first", "Perform depth first search,",
        "breadth-first", "Perform breadth first search,",
        "dynamic", "Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 "
        "integer feasible solutions have been found).",
        "best-guess", "choose node with smallest guessed integer solution",
        "Choose the strategy for selecting the next node to be processed.");
    roptions->setOptionExtraInfo("node_comparison", 63);

    roptions->AddStringOption5("tree_search_strategy",
        "Pick a strategy for traversing the tree",
        "probed-dive",
        "top-node"," Always pick the top node as sorted by the node comparison function",
        "dive","Dive in the tree if possible, otherwise pick top node as sorted by the tree comparison function.",
        "probed-dive","Dive in the tree exploring two children before continuing the dive at each level.",
        "dfs-dive","Dive in the tree if possible doing a depth first search. "
        "Backtrack on leaves or when a prescribed depth is attained or "
        "when estimate of best possible integer feasible solution in subtree "
        "is worst than cutoff. "
        "Once a prescribed limit of backtracks is attained pick top node "
        "as sorted by the tree comparison function",
        "dfs-dive-dynamic","Same as dfs-dive but once enough solution are found switch to best-bound and if too many nodes switch to depth-first.",
        "All strategies can be used in conjunction with any of the node comparison functions. "
        "Options which affect dfs-dive are max-backtracks-in-dive and max-dive-depth. "
        "The dfs-dive won't work in a non-convex problem where objective does not decrease down branches."
                              );
    roptions->setOptionExtraInfo("tree_search_strategy", 63);

    roptions->AddLowerBoundedIntegerOption("number_strong_branch",
        "Choose the maximum number of variables considered for strong branching.",
        0,20,
        "Set the number of variables on which to do strong branching.");
    roptions->setOptionExtraInfo("number_strong_branch", 127);

 
    roptions->AddLowerBoundedIntegerOption
    ("number_before_trust",
     "Set the number of branches on a variable before its pseudo costs are to be believed "
     "in dynamic strong branching.",
     0,8,
     "A value of 0 disables pseudo costs.");
    roptions->setOptionExtraInfo("number_before_trust", 127);

    roptions->AddStringOption2("nlp_failure_behavior",
        "Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not "
        "able to guarantee optimality within the specified tolerances).",
        "stop",
        "stop", "Stop when failure happens.",
        "fathom", "Continue when failure happens.",
        "If set to \"fathom\", the algorithm will fathom the node when Ipopt fails to find a solution to the nlp "
        "at that node within the specified tolerances. "
        "The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.");
    roptions->setOptionExtraInfo("nlp_failure_behavior", 8);

    roptions->AddStringOption2("sos_constraints",
        "Whether or not to activate SOS constraints.",
        "enable",
        "enable","",
        "disable","",
        "(only type 1 SOS are supported at the moment)");
    roptions->setOptionExtraInfo("sos_constraints", 63);

#ifdef BONMIN_CURVATURE_BRANCHING
    roptions->AddStringOption10("variable_selection",
#else
    roptions->AddStringOption9("variable_selection",
#endif
        "Chooses variable selection strategy",
        "strong-branching",
        "most-fractional", "Choose most fractional variable",
        "strong-branching", "Perform strong branching",
        "reliability-branching", "Use reliability branching",
#ifdef BONMIN_CURVATURE_BRANCHING
        "curvature-estimator", "Use curvature estimation to select branching variable",
#endif
        "qp-strong-branching", "Perform strong branching with QP approximation",
        "lp-strong-branching", "Perform strong branching with LP approximation",
        "nlp-strong-branching", "Perform strong branching with NLP approximation",
        "osi-simple", "Osi method to do simple branching",
        "osi-strong", "Osi method to do strong branching",
        "random", "Method to choose branching variable randomly");

    roptions->setOptionExtraInfo("variable_selection", 27);

    roptions->AddLowerBoundedIntegerOption("num_cut_passes",
        "Set the maximum number of cut passes at regular nodes of the branch-and-cut.",
        0,1,
        "");
    roptions->setOptionExtraInfo("num_cut_passes", 19);

    roptions->AddLowerBoundedIntegerOption("num_cut_passes_at_root",
        "Set the maximum number of cut passes at regular nodes of the branch-and-cut.",
        0,20,
        "");
    roptions->setOptionExtraInfo("num_cut_passes_at_root", 19);

    roptions->AddStringOption2("enable_dynamic_nlp",
               "Enable dynamic linear and quadratic rows addition in nlp",
               "no",
               "no", "",
               "yes", "",
               "");
    roptions->setOptionExtraInfo("enable_dynamic_nlp", 19);


    //roptions->SetRegisteringCategory("Debugging",RegisteredOptions::UndocumentedCategory);
    roptions->AddStringOption2("read_solution_file",
               "Read a file with the optimal solution to test if algorithms cuts it.",
               "no",
               "no","",
               "yes","",
               "For Debugging purposes only.");
    roptions->setOptionExtraInfo("enable_dynamic_nlp", 8);

    /* Branching options.*/
    LpBranchingSolver::registerOptions(roptions);

#ifdef COIN_HAS_FILTERSQP
    FilterSolver::registerOptions(roptions);
    BqpdSolver::registerOptions(roptions);
#endif
    CbcDiver::registerOptions(roptions);
    CbcDfsDiver::registerOptions(roptions);
    BonChooseVariable::registerOptions(roptions);
  }


  /** Initialize the options and the journalist.*/
  void
  BabSetupBase::initializeOptionsAndJournalist()
  {
    options_ = new Ipopt::OptionsList();

    journalist_= new Ipopt::Journalist();
    roptions_ = new Bonmin::RegisteredOptions();

    try {
      Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
        journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);

      options_->SetJournalist(journalist_);
      options_->SetRegisteredOptions(GetRawPtr(roptions_));
    }
    catch (Ipopt::IpoptException &E) {
      E.ReportException(*journalist_);
      throw E;
    }
    catch (std::bad_alloc) {
      journalist_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "\n Not enough memory .... EXIT\n");
      throw -1;
    }
#ifndef NO_CATCH_ALL
    catch (...) {
      Ipopt::IpoptException E("Uncaught exception in FilterSolver::FilterSolver()",
          "BonFilterSolver.cpp",-1);
      throw E;
    }
#endif

    registerOptions();
  }

  /** Get the options from given fileName */
  void
  BabSetupBase::readOptionsFile(std::string fileName)
  {
    if (GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
      initializeOptionsAndJournalist();
    std::ifstream is;
    if (fileName != "") {
      try {
        is.open(fileName.c_str());
      }
      catch (std::bad_alloc) {
        journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
        throw -1;
      }
#ifndef NO_CATCH_ALL
      catch (...) {
        Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
        E.ReportException(*journalist_);
        throw -1;
      }
#endif
    }
    readOptionsStream(is);
    if (is) {
      is.close();
    }
  }

  /** Get the options from long string containing all.*/
  void
  BabSetupBase::readOptionsString(std::string opt_string)
  {
    if (GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
      initializeOptionsAndJournalist();
    std::stringstream is(opt_string.c_str());
    readOptionsStream(is);
  }


  void
  BabSetupBase::readOptionsStream(std::istream& is)
  {
    if (GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
      initializeOptionsAndJournalist();
    if (is.good()) {
      try {
        options_->ReadFromStream(*journalist_, is);
      }
      catch (Ipopt::IpoptException &E) {
        E.ReportException(*journalist_);
        throw E;
      }
    }
    mayPrintDoc();
    readOptions_=true;
  }

  /** May print documentation of options if options print_options_documentation is set to yes.*/
  void
  BabSetupBase::mayPrintDoc()
  {
    bool print_options_documentation;
    options_->GetBoolValue("print_options_documentation",
        print_options_documentation, "");
    if (print_options_documentation) {
      std::list<std::string> categories;
      categories.push_back("Algorithm choice");
      categories.push_back("Branch-and-bound options");
      categories.push_back("ECP cuts generation");
      categories.push_back("Feasibility checker using OA cuts");
      categories.push_back("MILP Solver");
      categories.push_back("MILP cutting planes in hybrid algorithm");
      categories.push_back("Primal Heuristics");
      categories.push_back("NLP interface");
      categories.push_back("NLP solution robustness");
      categories.push_back("NLP solves in hybrid algorithm");
      categories.push_back("Nonconvex problems");
      categories.push_back("Outer Approximation Decomposition (B-OA)");
      categories.push_back("Outer Approximation cuts generation");
      categories.push_back("Output and Loglevel");
      categories.push_back("Strong branching setup");
      // Undocumented categories
      categories.push_back("Diving options");
      categories.push_back("ECP based strong branching");
      categories.push_back("Primal Heuristics (undocumented)");
      categories.push_back("Outer Approximation strengthening");
#ifdef COIN_HAS_FILTERSQP
      categories.push_back("FilterSQP options");
#endif
      //    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
      roptions_->OutputOptionDocumentation(*(journalist_),categories);
    }
  }

  void
  BabSetupBase::setPriorities()
  {
    const int * priorities = nonlinearSolver()->getPriorities();
    const double * upPsCosts = nonlinearSolver()->getUpPsCosts();
    const int * directions = nonlinearSolver()->getBranchingDirections();
    bool hasPseudo = (upPsCosts!=NULL);
    if (priorities == NULL && directions && NULL && hasPseudo)
      return;
    int n = nonlinearSolver()->numberObjects();
    OsiObject ** objects = nonlinearSolver()->objects();
    for (int i = 0 ; i < n; i++) {
      OsiObject2 * object = dynamic_cast<OsiObject2 *>(objects[i]);
      int iCol = objects[i]->columnNumber();
      if (iCol < 0) {
        throw CoinError("BabSetupBase","setPriorities",
            "Don't know how to set priority for non-column object.");
      }
      if (priorities) {
        objects[i]->setPriority(priorities[iCol]);
      }
      if (directions) {
        if (object == NULL) {
          throw CoinError("BabSetupBase","setPriorities",
              "Don't know how to set preferred way for object.");
        }
        object->setPreferredWay(directions[iCol]);
      }
      if (upPsCosts) {
        throw CoinError("BabSetupBase","setPriorities",
            "Can not handle user set pseudo-costs with OsiObjects\n"
            "You should use one of the Cbc branching rules:\n"
            "most-fractional or strong-branching.");
      }
    }
  }

  void
  BabSetupBase::addSos()
  {

    // pass user set Sos constraints (code inspired from CoinSolve.cpp)
    const TMINLP::SosInfo * sos = nonlinearSolver()->model()->sosConstraints();
    if (!getIntParameter(BabSetupBase::DisableSos) && sos && sos->num > 0) //we have some sos constraints
    {
      const int & numSos = sos->num;
      OsiObject ** objects = new OsiObject*[numSos];
      const int * starts = sos->starts;
      const int * indices = sos->indices;
      const char * types = sos->types;
      const double * weights = sos->weights;
      bool hasPriorities = false;
      const int * varPriorities = nonlinearSolver()->getPriorities();
      int numberObjects =  nonlinearSolver()->numberObjects();
      if (varPriorities)
      {
        for (int i = 0 ; i < numberObjects ; i++) {
          if (varPriorities[i]) {
            hasPriorities = true;
            break;
          }
        }
      }
      const int * sosPriorities = sos->priorities;
      if (sosPriorities)
      {
        for (int i = 0 ; i < numSos ; i++) {
          if (sosPriorities[i]) {
            hasPriorities = true;
            break;
          }
        }
      }
      for (int i = 0 ; i < numSos ; i++)
      {
        int start = starts[i];
        int length = starts[i + 1] - start;
        objects[i] = new OsiSOS(nonlinearSolver(), length, &indices[start],
            &weights[start], (int) types[i]);

        objects[i]->setPriority(10);
        if (hasPriorities && sosPriorities && sosPriorities[i]) {
          objects[i]->setPriority(sosPriorities[i]);
        }
      }
      nonlinearSolver()->addObjects(numSos, objects);
      for (int i = 0 ; i < numSos ; i++)
        delete objects[i];
      delete [] objects;
    }
  }


}/* End namespace Bonmin.*/

