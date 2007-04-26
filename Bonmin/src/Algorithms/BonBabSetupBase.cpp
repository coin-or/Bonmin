// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007


#include "BonBabSetupBase.hpp"
#include <fstream>
#include <sstream>


namespace Bonmin{
  int BabSetupBase::defaultIntParam_[BabSetupBase::NumberIntParam] = {
  1 /* BabLogLevel*/,
  100 /* BabLogInterval*/,
  2 /* MaxFailures.*/,
  0 /* FailureBehavior.*/,
  0 /* MaxInfeasible*/,
  5 /*NumberStrong*/,
  2 /* MinReliability*/,
  1 /*NumEcpRoundsStrong*/,
  INT_MAX /* MaxNodes*/,
  INT_MAX /* MaxSolutions*/,
  INT_MAX /* MaxIterations*/,
  0 /* SpecialOption*/,
  0 /* DisableSos.*/,
};


double BabSetupBase::defaultDoubleParam_[BabSetupBase::NumberDoubleParam] = {
  0 /* CutoffDecr*/,
  DBL_MAX /* Cutoff*/,
  0 /* AllowableGap*/,
  0 /*AllowableFractionGap*/,
  1e-09 /*IntTol*/,
  DBL_MAX /* MaxTime*/,
};

  BabSetupBase::BabSetupBase():
  nonlinearSolver_(NULL),
  continuousSolver_(NULL),
  cutGenerators_(),
  heuristics_(),
  branchingMethod_(NULL),
  nodeSelectionMethod_(),
  journalist_(NULL),
  options_(NULL),
roptions_(NULL),
readOptions_(false)
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
nodeSelectionMethod_(other.nodeSelectionMethod_),
journalist_(other.journalist_),
options_(NULL),
roptions_(other.roptions_),
readOptions_(other.readOptions_)
{
  if(other.nonlinearSolver_){
    nonlinearSolver_ = (OsiTMINLPInterface *) other.nonlinearSolver_->clone();}
  if(other.continuousSolver_){
    continuousSolver_ = other.continuousSolver_->clone();}
 
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
nodeSelectionMethod_(),
readOptions_(false)
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
nodeSelectionMethod_(),
journalist_(NULL),
options_(NULL),
roptions_(NULL),
readOptions_(false)
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
nodeSelectionMethod_(),
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
    delete nonlinearSolver_;
  delete continuousSolver_;
  delete branchingMethod_;
  for(CuttingMethods::iterator i = cutGenerators_.begin() ; i != cutGenerators_.end() ; i++){
    delete i->cgl;
    i->cgl = NULL;
  }
  
  for(HeuristicMethods::iterator i = heuristics_.begin() ; i != heuristics_.end() ; i++){
    delete *i;
  }
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
  options->GetIntegerValue("number_ecp_rounds_strong",intParam_[NumEcpRoundsStrong],"bonmin.");
  options->GetIntegerValue("node_limit",intParam_[MaxNodes],"bonmin.");
  options->GetIntegerValue("solution_limit",intParam_[MaxSolutions],"bonmin.");
  options->GetIntegerValue("iteration_limit",intParam_[MaxIterations],"bonmin.");
  options->GetEnumValue("sos_constraints",intParam_[DisableSos],"bonmin.");
  
  options->GetNumericValue("cutoff_decr",doubleParam_[CutoffDecr],"bonmin.");
  options->GetNumericValue("cutoff",doubleParam_[Cutoff],"bonmin.");
  options->GetNumericValue("allowable_gap",doubleParam_[AllowableGap],"bonmin.");
  options->GetNumericValue("allowable_fraction_gap",doubleParam_[AllowableFractionGap],"bonmin.");
  options->GetNumericValue("integer_tolerance",doubleParam_[IntTol],"bonmin.");
  options->GetNumericValue("time_limit", doubleParam_[MaxTime],"bonmin.");

  int ival;
  options->GetEnumValue("nodeselect_stra",ival,"bonmin.");
  nodeSelectionMethod_ = NodeSelectionStrategy(ival);
  
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
BabSetupBase::registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
  OsiTMINLPInterface::registerOptions(roptions);
  /* BabSetup options.*/
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
  
  roptions->SetRegisteringCategory("bonmin branch-and-bound options");
  
  roptions->AddLowerBoundedNumberOption("time_limit",
                                        "Set the global maximum computation time (in secs) for the algorithm.",
                                        0.,1,1e10,
                                        "");
  
  roptions->AddLowerBoundedIntegerOption("node_limit",
                                         "Set the maximum number of nodes explored in the branch-and-bound search.",
                                         0,INT_MAX,
                                         "");
  
  roptions->AddLowerBoundedIntegerOption("iteration_limit",
                                         "Set the cummulated maximum number of iteration in the algorithm used to process nodes continuous relaxations in the branch-and-bound.",
                                         0,INT_MAX,
                                         "value 0 deactivates option.");
  
  roptions->AddLowerBoundedIntegerOption("solution_limit",
                                         "Abort after that much integer feasible solution have been found by algorithm",
                                         0,INT_MAX,
                                         "value 0 deactivates option");
  
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
                             "dynamic",
                             "best-bound", "choose node whith the smallest bound,",
                             "depth-first", "Perform depth first search,",
                             "breadth-first", "Perform breadth first search,",
                             "dynamic", "Cbc dynamic strategy (starts with a depth first search and turn to best bound after 3 "
                             "integer feasible solutions have been found).",
                             "Choose the strategy for selecting the next node to be processed.");
  
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
  
  roptions->AddStringOption2("nlp_failure_behavior",
                             "Set the behavior when an NLP or a series of NLP are unsolved by Ipopt (we call unsolved an NLP for which Ipopt is not "
                             "able to guarantee optimality within the specified tolerances).",
                             "stop",
                             "stop", "Stop when failure happens.",
                             "fathom", "Continue when failure happens.",
                             "If set to \"fathom\", the algorithm will fathom the node when Ipopt fails to find a solution to the nlp "
                             "at that node whithin the specified tolerances. "
                             "The algorithm then becomes a heuristic, and the user will be warned that the solution might not be optimal.");
  
}


/** Initialize the options and the journalist.*/
void 
BabSetupBase::initializeOptionsAndJournalist(){
  options_ = new Ipopt::OptionsList();
  
  journalist_= new Ipopt::Journalist();
  roptions_ = new Ipopt::RegisteredOptions();
  
  try{
    Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
    journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
    
    options_->SetJournalist(journalist_);
    options_->SetRegisteredOptions(roptions_);
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
  std::cout<<"Reading options...."<<std::endl;
  if(is.good()){
    options_->ReadFromStream(*journalist_, is);
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
    //    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
    roptions_->OutputOptionDocumentation(*(journalist_),categories);
  }
}


}/* End namespace Bonmin.*/

