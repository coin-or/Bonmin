// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/19/2007


#include "BonCbc.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonAuxInfos.hpp"
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"

#include "CbcCompareUser.hpp"
#include "CbcCompareActual.hpp"

#include "CbcBranchUser.hpp"

// sets cutoff a bit above real one, to avoid single-point feasible sets
#define CUTOFF_TOL 1e-6

// Code to enable user interuption
static CbcModel * currentBranchModel = NULL; //pointer to the main b&b
Bonmin::OACutGenerator2 * currentOA = NULL; //pointer to the OA generator
CbcModel * OAModel; // pointer to the submip if using Cbc
bool BonminAbortAll;

#define SIGNAL
#ifdef SIGNAL
#include "CoinSignal.hpp"

extern "C"
{
  
  static bool BonminInteruptedOnce =false;
  static void signal_handler(int whichSignal) {
    if(BonminInteruptedOnce)
    {
      std::cerr<<"User forced interuption"<<std::endl;
      exit(0);
    }
    if (currentBranchModel!=NULL)
      currentBranchModel->setMaximumNodes(0); // stop at next node
    if (OAModel!=NULL)
      OAModel->setMaximumNodes(0); // stop at next node
    if (currentOA!=NULL)
      currentOA->parameter().maxLocalSearchTime_ = 0.; // stop OA
    BonminAbortAll = true;
    BonminInteruptedOnce = true;
    return;
  }
}
#endif

namespace Bonmin
{
  int usingCouenne = 0;
  
  /** Constructor.*/
  Bab2::Bab2():
  bestSolution_(NULL),
  mipStatus_(),
  bestObj_(1e200),
  bestBound_(-1e200),
  continuousRelaxation_(-COIN_DBL_MAX),
  numNodes_(0),
  mipIterationCount_(0),
  model_()
  {}
  
  /** Destructor.*/
  Bab2::~Bab2()
  {
    if (bestSolution_) delete [] bestSolution_;
    bestSolution_ = NULL;
  }
  
  /**operator() performs the branchAndBound*/
  void 
  Bab2::operator()(BabSetupBase & s)
  {
    branchAndBound(s);
  }
  
  /** Perform a branch-and-bound on given setup.*/
  void
  Bab2::branchAndBound(BabSetupBase & s)
  {
    /* Put a link to this into solver.*/
    OsiBabSolver *  babInfo = dynamic_cast<OsiBabSolver *>(s.continuousSolver()->getAuxiliaryInfo());
    Bonmin::BabInfo bonBabInfo(*babInfo);
    bonBabInfo.setBabPtr(this);
    s.continuousSolver()->setAuxiliaryInfo(&bonBabInfo);
    OsiSolverInterface * solver = s.continuousSolver()->clone();
    model_.passInMessageHandler(s.continuousSolver()->messageHandler());
    model_.assignSolver(solver, true);

  //  s.continuousSolver() = model_.solver();
    
    
    if(s.continuousSolver()->objects()!=NULL){
      model_.addObjects(s.continuousSolver()->numberObjects(),s.continuousSolver()->objects());
    }
    
    
    
    int specOpt = s.getIntParameter(BabSetupBase::SpecialOption);
    if(specOpt){
      model_.setSpecialOptions(specOpt);
      if(specOpt==16){
        CbcNlpStrategy strat(s.getIntParameter(BabSetupBase::MaxFailures), s.getIntParameter(BabSetupBase::MaxInfeasible), s.getIntParameter(BabSetupBase::FailureBehavior));
        model_.setStrategy(strat);
      }
    }
    
  
  
    model_.setMaximumCutPasses(s.getIntParameter(BabSetupBase::NumCutPasses));
   
    //Setup cutting plane methods
    for(BabSetupBase::CuttingMethods::iterator i = s.cutGenerators().begin() ;
        i != s.cutGenerators().end() ; i++){

      OaDecompositionBase * oa = dynamic_cast<OaDecompositionBase *>(i->cgl);
      if (oa && oa->reassignLpsolver())
        oa->assignLpInterface(model_.solver());
      if(i->atSolution)
        model_.addCutGenerator(i->cgl,i->frequency,i->id.c_str(), false, true);
      else
        model_.addCutGenerator(i->cgl,i->frequency,i->id.c_str());
    }
    
    for(BabSetupBase::HeuristicMethods::iterator i = s.heuristics().begin() ; 
        i != s.heuristics().end() ; i++){
      CbcHeuristic * heu = *i;
      heu->setModel(&model_);
      model_.addHeuristic(*i);
    }
    
    
    //need to record solver logLevel here
    int logLevel = s.continuousSolver()->messageHandler()->logLevel();
    
    //Set true branch-and-bound parameters
    model_.setLogLevel(s.getIntParameter(BabSetupBase::BabLogLevel));    
    
    // Put back solver logLevel
    model_.solver()->messageHandler()->setLogLevel(logLevel);
    
    model_.setPrintFrequency(s.getIntParameter(BabSetupBase::BabLogInterval));
        
    bool ChangedObject = false;
    //Pass over user set branching priorities to Cbc
    if(s.continuousSolver()->objects()==NULL)
    {
      const OsiTMINLPInterface * nlpSolver = s.nonlinearSolver();
      //set priorities, prefered directions...
      const int * priorities = nlpSolver->getPriorities();
      const double * upPsCosts = nlpSolver->getUpPsCosts();
      const double * downPsCosts = nlpSolver->getDownPsCosts();
      const int * directions = nlpSolver->getBranchingDirections();
      bool hasPseudo = (upPsCosts!=NULL);
      model_.findIntegers(true,hasPseudo);
      OsiObject ** simpleIntegerObjects = model_.objects();
      int numberObjects = model_.numberObjects();
      if(priorities != NULL || directions != NULL || hasPseudo)
      {
        ChangedObject = true;
        for (int i = 0 ; i < numberObjects ; i++)
        {
          CbcObject * object = dynamic_cast<CbcObject *>
          (simpleIntegerObjects[i]);
          int iCol = object->columnNumber();
          if (priorities)
            object->setPriority(priorities[iCol]);
          if (directions)
            object->setPreferredWay(directions[iCol]);
          if (upPsCosts) {
            CbcSimpleIntegerPseudoCost * pscObject =
            dynamic_cast<CbcSimpleIntegerPseudoCost*> (object);
            pscObject->setUpPseudoCost(upPsCosts[iCol]);
            pscObject->setDownPseudoCost(downPsCosts[iCol]);
          }
        }
      }
    }
    
    // Now pass user set Sos constraints (code inspired from CoinSolve.cpp)
    const TMINLP::SosInfo * sos = s.nonlinearSolver()->model()->sosConstraints();
    if (!s.getIntParameter(BabSetupBase::DisableSos) && sos && sos->num > 0) //we have some sos constraints
    {
      ChangedObject = true;
      const OsiTMINLPInterface * nlpSolver = s.nonlinearSolver();
      const int & numSos = sos->num;
      CbcObject ** objects = new CbcObject*[numSos];
      const int * starts = sos->starts;
      const int * indices = sos->indices;
      const char * types = sos->types;
      const double * weights = sos->weights;
      //verify if model has user set priorities
      bool hasPriorities = false;
      const int * varPriorities = nlpSolver->getPriorities();
      int numberObjects = model_.numberObjects();
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
        objects[i] = new CbcSOS(&model_, length, &indices[start],
                                &weights[start], i, types[i]);
        
        objects[i]->setPriority(10);
        if (hasPriorities && sosPriorities && sosPriorities[i]) {
          objects[i]->setPriority(sosPriorities[i]);
        }
      }
      model_.addObjects(numSos, objects);
      for (int i = 0 ; i < numSos ; i++)
        delete objects[i];
      delete [] objects;
    }
    
    replaceIntegers(model_.objects(), model_.numberObjects());
    
    
    model_.setDblParam(CbcModel::CbcCutoffIncrement, s.getDoubleParameter(BabSetupBase::CutoffDecr));
    
    model_.setCutoff(s.getDoubleParameter(BabSetupBase::Cutoff) + CUTOFF_TOL);
    
    model_.setDblParam(CbcModel::CbcAllowableGap, s.getDoubleParameter(BabSetupBase::AllowableGap));
    model_.setDblParam(CbcModel::CbcAllowableFractionGap, s.getDoubleParameter(BabSetupBase::AllowableFractionGap));
    
    // Definition of node selection strategy

    if (s.nodeSelectionMethod()==BabSetupBase::bestBound) {
      CbcCompareObjective compare;
      model_.setNodeComparison(compare);
    }
    else if (s.nodeSelectionMethod()==BabSetupBase::DFS) {
      CbcCompareDepth compare;
      model_.setNodeComparison(compare);
    }
    else if (s.nodeSelectionMethod()==BabSetupBase::BFS) {
      CbcCompareUser compare;
      compare.setWeight(0.0);
      model_.setNodeComparison(compare);
    }
    else if (s.nodeSelectionMethod()==BabSetupBase::dynamic) {
      CbcCompareUser compare;
      model_.setNodeComparison(compare);
    }
    
    model_.setNumberStrong(s.getIntParameter(BabSetupBase::NumberStrong));
    
    model_.setNumberBeforeTrust(s.getIntParameter(BabSetupBase::MinReliability));
    
    model_.setNumberPenalties(8);
    
    model_.setDblParam(CbcModel::CbcMaximumSeconds, s.getDoubleParameter(BabSetupBase::MaxTime));
    
    model_.setMaximumNodes(s.getIntParameter(BabSetupBase::MaxNodes));
    
    model_.setMaximumSolutions(s.getIntParameter(BabSetupBase::MaxSolutions));
    
    model_.setIntegerTolerance(s.getDoubleParameter(BabSetupBase::IntTol));
    
    
    // Redundant definition of default branching (as Default == User)
    CbcBranchUserDecision branch;
    if(s.branchingMethod() != NULL){
      s.branchingMethod()->setSolver(model_.solver());
      s.branchingMethod()->setNumberStrong(model_.numberStrong());
      OsiChooseStrong * strong = dynamic_cast<OsiChooseStrong *>(s.branchingMethod());
      if(strong)
        strong->setNumberBeforeTrusted(model_.numberBeforeTrust());
      branch.setChooseMethod(*s.branchingMethod());
    }
    
    
    model_.setBranchingMethod(&branch);
    
    //Get the time and start.
    model_.initialSolve();
    
    continuousRelaxation_ =model_.solver()->getObjValue();
    if(specOpt==16)//Set warm start point for Ipopt
    {
      const double * colsol = model_.solver()->getColSolution();
      const double * duals = model_.solver()->getRowPrice();
      model_.solver()->setColSolution(colsol);
      model_.solver()->setRowPrice(duals);
    }
    
#ifdef SIGNAL
    CoinSighandler_t saveSignal=SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
#endif
    
    currentBranchModel = &model_;
    model_.branchAndBound();
    
    
    numNodes_ = model_.getNodeCount();
    bestObj_ = model_.getObjValue();
    bestBound_ = model_.getBestPossibleObjValue();
    mipIterationCount_ = model_.getIterationCount();
    
    bool hasFailed = false;
    if(specOpt==16)//Did we continue branching on a failure
    {
      CbcNlpStrategy * nlpStrategy = dynamic_cast<CbcNlpStrategy *>(model_.strategy());
      if (nlpStrategy)
        hasFailed = nlpStrategy->hasFailed();
      else
        throw -1;
    }
    else
      hasFailed = s.nonlinearSolver()->hasContinuedOnAFailure();
    
    
    if (hasFailed) {
      std::cout<<"************************************************************"<<std::endl
      <<"WARNING : Optimization failed on an NLP during optimization"
      <<"\n (no optimal value found within tolerances)."<<std::endl
      <<"Optimization was not stopped because option \n"
      <<"\"nlp_failure_behavior\" has been set to fathom but"
      <<" beware that reported solution may not be optimal"<<std::endl
      <<"************************************************************"<<std::endl;
    }
    TMINLP::SolverReturn status;
    
    if(model_.numberObjects()==0){
      if(bestSolution_)
        delete [] bestSolution_;
      bestSolution_ = new double[s.nonlinearSolver()->getNumCols()];
      CoinCopyN(s.nonlinearSolver()->getColSolution(), s.nonlinearSolver()->getNumCols(),
                bestSolution_);
      bestObj_ = bestBound_ = s.nonlinearSolver()->getObjValue();
    }
    
    if (model_.bestSolution()) {
      if (bestSolution_)
        delete [] bestSolution_;
      bestSolution_ = new double[s.nonlinearSolver()->getNumCols()];
      CoinCopyN(model_.bestSolution(), s.nonlinearSolver()->getNumCols(), bestSolution_);
    }
    if (model_.status() == 0) {
      if (bestSolution_){
        status = TMINLP::SUCCESS;
        mipStatus_ = FeasibleOptimal;
      }
      else {
        status = TMINLP::INFEASIBLE;
        mipStatus_ = ProvenInfeasible;
      }
    }
    else if(model_.status() == 1){
      status = TMINLP::LIMIT_EXCEEDED;
      if (bestSolution_){
        mipStatus_ = Feasible;
      }
      else {
        mipStatus_ = NoSolutionKnown;
      }
    }
    else if(model_.status()==2){
      status = TMINLP::MINLP_ERROR;
    }
    s.nonlinearSolver()->model()->finalize_solution(status, s.nonlinearSolver()->getNumCols(), bestSolution_,
                                          bestObj_);
  }
  
  
  /** return the best known lower bound on the objective value*/
  double
  Bab2::bestBound()
  {
    if (mipStatus_ == FeasibleOptimal) return bestObj_;
    else if (mipStatus_ == ProvenInfeasible) return 1e200;
    else return bestBound_;
  }
}
