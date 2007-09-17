// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
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

namespace Bonmin{
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
  20 /* numCutPassesAtRoot.*/
};


double BabSetupBase::defaultDoubleParam_[BabSetupBase::NumberDoubleParam] = {
  0 /* CutoffDecr*/,
  COIN_DBL_MAX /* Cutoff*/,
  0 /* AllowableGap*/,
  0 /*AllowableFractionGap*/,
  1e-09 /*IntTol*/,
  COIN_DBL_MAX /* MaxTime*/,
};

  BabSetupBase::BabSetupBase():
  nonlinearSolver_(NULL),
  continuousSolver_(NULL),
  cutGenerators_(),
  heuristics_(),
  branchingMethod_(NULL),
  nodeComparisonMethod_(),
  treeTraversalMethod_(),
  journalist_(NULL),
  options_(NULL),
  roptions_(NULL),
  readOptions_(false),
  lpMessageHandler_(NULL)
{
  CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
  CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
}

/** Copy constructor. */
BabSetupBase::BabSetupBase(const BabSetupBase & other):
nonlinearSolver_(NULL),
continuousSolver_(NULL),
cutGenerators_(),
heuristics_(),
branchingMethod_(),
nodeComparisonMethod_(other.nodeComparisonMethod_),
treeTraversalMethod_(other.treeTraversalMethod_),
journalist_(other.journalist_),
options_(NULL),
roptions_(other.roptions_),
readOptions_(other.readOptions_)
{
  if(other.nonlinearSolver_){
    nonlinearSolver_ = static_cast<OsiTMINLPInterface *>(other.nonlinearSolver_->clone());}
  if(other.continuousSolver_){
    continuousSolver_ = other.continuousSolver_->clone();}
  if(other.lpMessageHandler_){
    lpMessageHandler_ = other.lpMessageHandler_->clone();}
  continuousSolver_->passInMessageHandler(lpMessageHandler_);
  for(CuttingMethods::const_iterator i = other.cutGenerators_.begin() ; i != other.cutGenerators_.end() ; i++){
    cutGenerators_.push_back(*i);
    cutGenerators_.back().cgl = cutGenerators_.back().cgl->clone();
  }
  
  for(HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++){
    heuristics_.push_back((*i)->clone());
  }
  branchingMethod_ = other.branchingMethod_->clone();
  if(IsValid(other.options_)){
    options_ = new OptionsList;
    *options_ = *other.options_;
  }
  CoinCopyN(other.intParam_, NumberIntParam, intParam_);
  CoinCopyN(other.doubleParam_, NumberDoubleParam, doubleParam_);
}

BabSetupBase::BabSetupBase(Ipopt::SmartPtr<TMINLP> tminlp):
nonlinearSolver_(NULL),
continuousSolver_(NULL),
cutGenerators_(),
heuristics_(),
branchingMethod_(NULL),
nodeComparisonMethod_(),
treeTraversalMethod_(),
readOptions_(false),
lpMessageHandler_(NULL)
{ 
  CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
  CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
  use(tminlp);
}




void 
BabSetupBase::use(Ipopt::SmartPtr<TMINLP> tminlp){
  readOptionsFile();
  assert(IsValid(tminlp));
  nonlinearSolver_ = new OsiTMINLPInterface;
  nonlinearSolver_->initialize(roptions_, options_, journalist_, tminlp);
}

BabSetupBase::BabSetupBase(const OsiTMINLPInterface& nlp):
nonlinearSolver_(NULL),
continuousSolver_(NULL),
cutGenerators_(),
heuristics_(),
branchingMethod_(NULL),
nodeComparisonMethod_(),
treeTraversalMethod_(),
journalist_(NULL),
options_(NULL),
roptions_(NULL),
readOptions_(false),
lpMessageHandler_(NULL)
{
  CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
  CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
  use(nlp);
}

void 
BabSetupBase::use(const OsiTMINLPInterface& nlp){
  nonlinearSolver_ = dynamic_cast<OsiTMINLPInterface *>(nlp.clone());
  options_ = nonlinearSolver_->options();
  roptions_ = nonlinearSolver_->regOptions();
  journalist_ = nonlinearSolver_->solver()->Jnlst();
  readOptions_ = true;
}

BabSetupBase::BabSetupBase( Ipopt::SmartPtr<TNLPSolver> app):
nonlinearSolver_(NULL),
continuousSolver_(NULL),
cutGenerators_(),
heuristics_(),
branchingMethod_(NULL),
nodeComparisonMethod_(),
treeTraversalMethod_(),
journalist_(app->Jnlst()),
options_(app->Options()),
roptions_(app->RegOptions()),
readOptions_(true)
{
  CoinCopyN(defaultIntParam_, NumberIntParam, intParam_);
  CoinCopyN(defaultDoubleParam_, NumberDoubleParam, doubleParam_);
}


BabSetupBase::~BabSetupBase(){
  if(nonlinearSolver_ != continuousSolver_)
  {
    delete nonlinearSolver_;
  }
  delete continuousSolver_;
  delete branchingMethod_;
  for(CuttingMethods::iterator i = cutGenerators_.begin() ; i != cutGenerators_.end() ; i++){
    delete i->cgl;
    i->cgl = NULL;
  }
  
  for(HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++){
    delete *i;
  }
  delete lpMessageHandler_;
}


void 
BabSetupBase::gatherParametersValues(Ipopt::SmartPtr<OptionsList> options){
  
  options->GetIntegerValue("bb_log_level",intParam_[BabLogLevel],"bonmin.");
  options->GetIntegerValue("bb_log_interval",intParam_[BabLogInterval],"bonmin.");
  options->GetIntegerValue("max_consecutive_failures",intParam_[MaxFailures],"bonmin.");
  options->GetEnumValue("nlp_failure_behavior",intParam_[FailureBehavior],"bonmin.");
  options->GetIntegerValue("max_consecutive_infeasible",intParam_[MaxInfeasible],"bonmin.");
  options->GetIntegerValue("number_strong_branch",intParam_[NumberStrong],"bonmin.");
  options->GetIntegerValue("number_before_trust",intParam_[MinReliability],"bonmin.");
  options->GetIntegerValue("node_limit",intParam_[MaxNodes],"bonmin.");
  options->GetIntegerValue("solution_limit",intParam_[MaxSolutions],"bonmin.");
  options->GetIntegerValue("iteration_limit",intParam_[MaxIterations],"bonmin.");
  options->GetEnumValue("sos_constraints",intParam_[DisableSos],"bonmin.");
  options->GetIntegerValue("num_cut_passes",intParam_[NumCutPasses],"bonmin.");
  options->GetIntegerValue("num_cut_passes_at_root",intParam_[NumCutPassesAtRoot],"bonmin.");
  
  options->GetNumericValue("cutoff_decr",doubleParam_[CutoffDecr],"bonmin.");
  options->GetNumericValue("cutoff",doubleParam_[Cutoff],"bonmin.");
  options->GetNumericValue("allowable_gap",doubleParam_[AllowableGap],"bonmin.");
  options->GetNumericValue("allowable_fraction_gap",doubleParam_[AllowableFractionGap],"bonmin.");
  options->GetNumericValue("integer_tolerance",doubleParam_[IntTol],"bonmin.");
  options->GetNumericValue("time_limit", doubleParam_[MaxTime],"bonmin.");

  int ival;
  options->GetEnumValue("node_comparison",ival,"bonmin.");
  nodeComparisonMethod_ = NodeComparison(ival);
  
  options->GetEnumValue("tree_search_strategy", ival, "bonmin.");
  treeTraversalMethod_ = TreeTraversal(ival);
  
  int varSelection;
  options->GetEnumValue("varselect_stra",varSelection,"bonmin.");
  // Set branching strategy
  if(varSelection == OsiTMINLPInterface::MOST_FRACTIONAL){
    intParam_[MinReliability] = 0;
    intParam_[NumberStrong] = 0;
  }
  else if(varSelection == OsiTMINLPInterface::STRONG_BRANCHING){
    intParam_[MinReliability] = 0;
  }
  else if(varSelection == OsiTMINLPInterface::RELIABILITY_BRANCHING){
    intParam_[MinReliability] = 10;
  }
}

void 
BabSetupBase::registerOptions(){
  registerAllOptions(roptions_);
}

void 
BabSetupBase::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
  OsiTMINLPInterface::registerOptions(roptions);
  /* BabSetup options.*/
  roptions->SetRegisteringCategory("bonmin output options", RegisteredOptions::BonminCategory);
  
  roptions->AddBoundedIntegerOption("bb_log_level",
                                    "specify main branch-and-bound log level.",
                                    0,5,1,
                                    "Set the level of output of the branch-and-bound : "
                                    "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                    );
  roptions->setOptionExtraInfo("bb_log_level", 11);
  
  roptions->AddLowerBoundedIntegerOption("bb_log_interval",
                                         "Interval at which node level output is printed.",
                                         0,100,
                                         "Set the interval (in terms of number of nodes) at which "
                                         "a log on node resolutions (consisting of lower and upper bounds) is given.");
  roptions->setOptionExtraInfo("bb_log_interval", 11);
  
  roptions->AddBoundedIntegerOption("lp_log_level",
                                    "specify LP log level.",
                                    0,4,0,
                                    "Set the level of output of the linear programming sub-solver in B-Hyb or B-QG : "
                                    "0 - none, 1 - minimal, 2 - normal low, 3 - normal high, 4 - verbose"
                                    );
  roptions->setOptionExtraInfo("lp_log_level", 3);
  
  roptions->SetRegisteringCategory("bonmin branch-and-bound options", RegisteredOptions::BonminCategory);
  
  roptions->AddLowerBoundedNumberOption("time_limit",
                                        "Set the global maximum computation time (in secs) for the algorithm.",
                                        0.,1,1e10,
                                        "");
  roptions->setOptionExtraInfo("time_limit", 31);
  
  roptions->AddLowerBoundedIntegerOption("node_limit",
                                         "Set the maximum number of nodes explored in the branch-and-bound search.",
                                         0,COIN_INT_MAX,
                                         "");
  roptions->setOptionExtraInfo("node_limit", 31);
  
  roptions->AddLowerBoundedIntegerOption("iteration_limit",
                                         "Set the cummulated maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound.",
                                         0,COIN_INT_MAX,
                                         "value 0 deactivates option.");
  roptions->setOptionExtraInfo("iteration_limit", 31);
  
  roptions->AddLowerBoundedIntegerOption("solution_limit",
                                         "Abort after that much integer feasible solution have been found by algorithm",
                                         0,COIN_INT_MAX,
                                         "value 0 deactivates option");
  roptions->setOptionExtraInfo("solution_limit", 31);
 
  roptions->AddBoundedNumberOption("integer_tolerance",
                                   "Set integer tolerance.",
                                   0.,1,.5,1,1e-06,
                                   "Any number within that value of an integer is considered integer.");
  roptions->setOptionExtraInfo("integer_tolerance", 31);
  
  roptions->AddBoundedNumberOption("allowable_gap",
                                   "Specify the value of absolute gap under which the algorithm stops.",
                                   -1.e20,0,1.e20,0,0.,
                                   "Stop the tree search when the gap between the objective value of the best known solution"
                                   " and the best bound on the objective of any solution is less than this.");
  roptions->setOptionExtraInfo("allowable_gap", 31);
  
  roptions->AddBoundedNumberOption("allowable_fraction_gap",
                                   "Specify the value of relative gap under which the algorithm stops.",
                                   -1.e20,0,1.e20,0,0.0,
                                   "Stop the tree search when the gap between the objective value of the best known solution "
                                   "and the best bound on the objective of any solution is less than this "
                                   "fraction of the absolute value of the best known solution value.");
  roptions->setOptionExtraInfo("allowable_fraction_gap", 31);
  
  roptions->AddBoundedNumberOption("cutoff",
                                   "Specify cutoff value.",
                                   -1e100,0,1e100,0,1e100,
                                   "cutoff should be the value of a feasible solution known by the user "
                                   "(if any). The algorithm will only look for solutions better than cutoof.");
  roptions->setOptionExtraInfo("cutoff", 31);
  
  
  roptions->AddBoundedNumberOption("cutoff_decr",
                                   "Specify cutoff decrement.",
                                   -1e10,0,1e10,0,1e-05,
                                   "Specify the amount by which cutoff is decremented below "
                                   "a new best upper-bound"
                                   " (usually a small postive value but in non-convex problems it may be a negative value).");
  roptions->setOptionExtraInfo("cutoff_decr", 31);
  
  
  roptions->AddStringOption5("node_comparison",
                             "Choose the node selection strategy.",
                             "dynamic",
                             "best-bound", "choose node whith the smallest bound,",
                             "depth-first", "Perform depth first search,",
                             "breadth-first", "Perform breadth first search,",
                             "dynamic", "Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 "
                             "integer feasible solutions have been found).",
			     "best-guess", "choose node with smallest guessed integer solution",
                             "Choose the strategy for selecting the next node to be processed.");
  roptions->setOptionExtraInfo("node_comparison", 31);

  roptions->AddStringOption4("tree_search_strategy",
			     "Pick a strategy for traversing the tree",
			     "top-node",
			     "top-node"," Always pick the top node as sorted by the node comparison function",
			     "dive","Dive in the tree if possible, otherwise pick top node as sorted by the tree comparison function",
			     "dfs-dive","Dive in the tree if possible doing a depth first search."
                                        "Backtrack on leaves or when a prescribed depth is attained or "
                                        "when estimate of best possible integer feasible solution in subtree "
                                        "is worst than cutoff. "
                                        "Once a prescribed limit of backtracks is attained pick top node "
                                        "as sorted by the tree comparison function",
			     "dfs-dive-dynamic","Same as dfs-dive but once enough solution are found switch to best-bound and if too many nodes switch to depth-first.",
			     "All strategies can be used in conjunction with any of the node comparison functions."
                             "Options which affect dfs-dive are max-backtracks-in-dive and max-dive-depth. "
                             "The dfs-dive won't work in a non-convex problem where objective does not decrease down branches."
                              ); 
  roptions->setOptionExtraInfo("tree_search_strategy", 31);
  
  roptions->AddLowerBoundedIntegerOption("number_strong_branch",
                                         "Choose the maximum number of variables considered for strong branching.",
                                         0,20,
                                         "Set the number of variables on which to do strong branching.");
  roptions->setOptionExtraInfo("number_strong_branch", 31);
  
  roptions->AddLowerBoundedIntegerOption
    ("number_before_trust",
     "Set the number of branches on a variable before its pseudo costs are to be believed "
     "in dynamic strong branching.",
     -1,8,
     "A value of -1 disables pseudo costs.");
  roptions->setOptionExtraInfo("number_before_trust", 31);
  
  roptions->AddStringOption2("nlp_failure_behavior",
                             "Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not "
                             "able to guarantee optimality within the specified tolerances).",
                             "stop",
                             "stop", "Stop when failure happens.",
                             "fathom", "Continue when failure happens.",
                             "If set to \"fathom\", the algorithm will fathom the node when Ipopt fails to find a solution to the nlp "
                             "at that node whithin the specified tolerances. "
                             "The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.");
  roptions->setOptionExtraInfo("nlp_failure_behavior", 31);
  
  roptions->AddStringOption2("sos_constraints",
                             "Wether or not to activate SOS constraints.",
                             "enable",
                             "enable","",
                             "disable","",
                             "(only type 1 SOS are supported at the moment)");
  roptions->setOptionExtraInfo("sos_constraints", 11);
  
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
  roptions->setOptionExtraInfo("varselect_stra", 15);
  
  roptions->AddLowerBoundedIntegerOption("num_cut_passes",
                                         "Set the maximum number of cut passes at regular nodes of the branch-and-cut.",
                                         0,1,
                                         "");
  roptions->setOptionExtraInfo("num_cut_passes", 7);

  roptions->AddLowerBoundedIntegerOption("num_cut_passes_at_root",
                                         "Set the maximum number of cut passes at regular nodes of the branch-and-cut.",
                                         0,20,
                                         "");
  roptions->setOptionExtraInfo("num_cut_passes_at_root", 7);

  roptions->SetRegisteringCategory("bonmin experimental options",  RegisteredOptions::BonminCategory);
  // Some options for the strong branching and pseudo-costs that we
  // still expored
  roptions->AddBoundedNumberOption("setup_pseudo_frac", "Proportion of strong branching list that has to be taken from most-integer-infeasible list.",
				   0., false, 1., false, 0.5);
  roptions->AddBoundedNumberOption("maxmin_crit_no_sol", "Weight towards minimum in of lower and upper branching estimates when no solution has been found yet.",
				   0., false, 1., false, 0.7);
  roptions->AddBoundedNumberOption("maxmin_crit_have_sol", "Weight towards minimum in of lower and upper branching estimates when a solution has been found.",
				   0., false, 1., false, 0.1);
  roptions->AddLowerBoundedIntegerOption("number_before_trust_list",
					 "Set the number of branches on a variable before its pseudo costs are to be believed during setup of strong branching candidate list.",
					 -1, 0, "The default value is that of \"number_before_trust\"");
  roptions->AddLowerBoundedIntegerOption("number_strong_branch_root",
					 "Maximum number of variables considered for strong branching in root node.",
					 0, COIN_INT_MAX, "");
  


    /* Branching options.*/
    LpBranchingSolver::registerOptions(roptions);

#ifdef COIN_HAS_FILTERSQP
    FilterSolver::registerOptions(roptions);
#endif
    CbcDiver::registerOptions(roptions);
    CbcDfsDiver::registerOptions(roptions);
}


/** Initialize the options and the journalist.*/
void 
BabSetupBase::initializeOptionsAndJournalist(){
  options_ = new Ipopt::OptionsList();
  
  journalist_= new Ipopt::Journalist();
  roptions_ = new Bonmin::RegisteredOptions();
  
  try{
    Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
    journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
    
    options_->SetJournalist(journalist_);
    options_->SetRegisteredOptions(GetRawPtr(roptions_));
  }
  catch (Ipopt::IpoptException &E){
    E.ReportException(*journalist_);
    throw E;
  }
  catch(std::bad_alloc){
    journalist_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "\n Not enough memory .... EXIT\n");
    throw -1;
  }
  catch(...){
    Ipopt::IpoptException E("Uncaught exception in FilterSolver::FilterSolver()",
                            "BonFilterSolver.cpp",-1);
    throw E;
  }
  
  registerOptions();
}

/** Get the options from given fileName */
void 
BabSetupBase::readOptionsFile(std::string fileName){
  if(GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
    initializeOptionsAndJournalist();
  std::ifstream is;
  if (fileName != "") {
    try {
      is.open(fileName.c_str());
    }
    catch(std::bad_alloc) {
      journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
      throw -1;
    }
    catch(...) {
      Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
      E.ReportException(*journalist_);
      throw -1;
    }
  }
  readOptionsStream(is);
  if (is) {
    is.close();
  }
}

/** Get the options from long string containing all.*/
void 
BabSetupBase::readOptionsString(std::string opt_string){
  if(GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
    initializeOptionsAndJournalist();
  std::stringstream is(opt_string.c_str());
  readOptionsStream(is);
}


void 
BabSetupBase::readOptionsStream(std::istream& is)
{
  if(GetRawPtr(options_) == NULL || GetRawPtr(roptions_) == NULL || GetRawPtr(journalist_) == NULL)
    initializeOptionsAndJournalist();
  if(is.good()){
    try{
      options_->ReadFromStream(*journalist_, is);
    }
    catch (Ipopt::IpoptException &E){
      E.ReportException(*journalist_);
      throw E;
    }
  }
  mayPrintDoc();
  readOptions_=true;
}

/** May print documentation of options if options print_options_documentation is set to yes.*/
void 
BabSetupBase::mayPrintDoc(){
  bool print_options_documentation;
  options_->GetBoolValue("print_options_documentation",
                         print_options_documentation, "");
  if (print_options_documentation) {
    std::list<std::string> categories;
    categories.push_back("bonmin branch-and-bound options");
    categories.push_back("bonmin options for robustness");
    categories.push_back("bonmin options for non-convex problems");
    categories.push_back("bonmin options : B-Hyb specific options");
#ifdef COIN_HAS_FILTERSQP
    categories.push_back("FilterSQP options");
#endif
    //    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
    roptions_->OutputOptionDocumentation(*(journalist_),categories);
  }
}


}/* End namespace Bonmin.*/

