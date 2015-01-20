// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006

#include <sstream>
#include <climits>

#include <algorithm>
#include "BonOaDecBase.hpp"


#include "BonminConfig.h"

#include "OsiClpSolverInterface.hpp"

#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "BonCbcLpStrategy.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#define CHECK_CPX_STAT(a,b) if(b) throw CoinError("Error in CPLEX call",__FILE__,a);

#endif
#include "BonCbc.hpp"
#include "BonSolverHelp.hpp"
//The following two are to interupt the solution of sub-mip through CTRL-C
extern CbcModel * OAModel;

namespace Bonmin {


  OaDecompositionBase::OaDecompositionBase(BabSetupBase &b, bool leaveSiUnchanged,
      bool reassignLpsolver):
      CglCutGenerator(),
      nlp_(b.nonlinearSolver()),
      s_(&b),
      lp_(NULL),
      objects_(NULL),
      nObjects_(0),
      nLocalSearch_(0),
      handler_(NULL),
      leaveSiUnchanged_(leaveSiUnchanged),
      reassignLpsolver_(reassignLpsolver),
      timeBegin_(0),
      numSols_(0),
      parameters_(),
      currentNodeNumber_(-1)
  {
    handler_ = new CoinMessageHandler();
    int logLevel;
    b.options()->GetIntegerValue("oa_log_level",logLevel,b.prefix());
    b.options()->GetNumericValue("oa_log_frequency",parameters_.logFrequency_,b.prefix());
    b.options()->GetNumericValue("allowable_fraction_gap", parameters_.gap_tol_, b.prefix());
    handler_ -> setLogLevel(logLevel);
    b.options()->GetIntegerValue("solution_limit", parameters_.maxSols_,b.prefix());

    messages_ = OaMessages();
    timeBegin_ = CoinCpuTime();
    b.options()->GetIntegerValue("milp_log_level",parameters_.subMilpLogLevel_,b.prefix());
    b.options()->GetNumericValue("cutoff_decr",parameters_.cbcCutoffIncrement_,b.prefix());
    b.options()->GetNumericValue("integer_tolerance",parameters_.cbcIntegerTolerance_,b.prefix());
    int ivalue;
    b.options()->GetEnumValue("add_only_violated_oa", ivalue,b.prefix());
    parameters_.addOnlyViolated_ = ivalue;
    b.options()->GetEnumValue("oa_cuts_scope", ivalue,b.prefix());
    parameters_.global_ = ivalue;
}

  OaDecompositionBase::OaDecompositionBase
  (const OaDecompositionBase & other)
      :
      CglCutGenerator(other),
      nlp_(other.nlp_),
      s_(other.s_),
      lp_(other.lp_),
      objects_(other.objects_),
      nObjects_(other.nObjects_),
      nLocalSearch_(0),
      messages_(other.messages_),
      leaveSiUnchanged_(other.leaveSiUnchanged_),
      reassignLpsolver_(other.reassignLpsolver_),
      timeBegin_(0),
      numSols_(other.numSols_),
      parameters_(other.parameters_),
      currentNodeNumber_(other.currentNodeNumber_)
  {
    timeBegin_ = CoinCpuTime();
    handler_ = other.handler_->clone();
  }
/// Constructor with default values for parameters
  OaDecompositionBase::Parameters::Parameters():
      global_(true),
      addOnlyViolated_(false),
      cbcCutoffIncrement_(1e-06),
      cbcIntegerTolerance_(1e-05),
      gap_tol_(1e-05),
      maxLocalSearch_(0),
      maxLocalSearchTime_(3600),
      subMilpLogLevel_(0),
      maxSols_(INT_MAX),
      logFrequency_(1000.),
      strategy_(NULL)
  {}

  /** Destructor.*/
  OaDecompositionBase::~OaDecompositionBase()
  {
    delete handler_;
  }


/// Constructor with default values for parameters
  OaDecompositionBase::Parameters::Parameters(const Parameters & other):
      global_(other.global_),
      addOnlyViolated_(other.addOnlyViolated_),
      cbcCutoffIncrement_(other.cbcCutoffIncrement_),
      cbcIntegerTolerance_(other.cbcIntegerTolerance_),
      gap_tol_(other.gap_tol_),
      maxLocalSearch_(other.maxLocalSearch_),
      maxLocalSearchTime_(other.maxLocalSearchTime_),
      subMilpLogLevel_(other.subMilpLogLevel_),
      maxSols_(other.maxSols_),
      logFrequency_(other.logFrequency_),
      strategy_(NULL)
  {
    if (other.strategy_)
      strategy_ = other.strategy_->clone();
  }



OaDecompositionBase::solverManip::solverManip
(OsiSolverInterface * si,
 bool saveNumRows,
 bool saveBasis,
 bool saveBounds,
 bool saveCutoff,
 bool resolve):
    si_(si),
    initialNumberRows_(-1),
    colLower_(NULL),
    colUpper_(NULL),
    warm_(NULL),
    cutoff_(DBL_MAX),
    deleteSolver_(false),
    objects_(NULL),
    nObjects_(0)
{
  getCached();
  if (saveNumRows)
    initialNumberRows_ = numrows_;
  if (saveBasis)
    warm_ = si->getWarmStart();
  if (saveBounds) {
    colLower_ = new double[numcols_];
    colUpper_ = new double[numcols_];
    CoinCopyN(si->getColLower(), numcols_ , colLower_);
    CoinCopyN(si->getColUpper(), numcols_ , colUpper_);
  }
  if (saveCutoff)
    si->getDblParam(OsiDualObjectiveLimit, cutoff_);
  si->messageHandler()->setLogLevel(0);
  if (resolve) si->resolve();
}


OaDecompositionBase::solverManip::solverManip
(const OsiSolverInterface & si):
    si_(NULL),
    initialNumberRows_(-1),
    colLower_(NULL),
    colUpper_(NULL),
    warm_(NULL),
    cutoff_(DBL_MAX),
    deleteSolver_(true),
    objects_(NULL),
    nObjects_(0)
{
  si_ = si.clone();
  getCached();
}

OaDecompositionBase::solverManip::~solverManip()
{
  if (warm_) delete warm_;
  if (colLower_) delete [] colLower_;
  if (colUpper_) delete [] colUpper_;
  if (deleteSolver_) delete si_;
}

void
OaDecompositionBase::solverManip::restore()
{
  if (initialNumberRows_ >= 0) {
    int nRowsToDelete = si_->getNumRows() - initialNumberRows_;
    int * rowsToDelete = new int[nRowsToDelete];
    for (int i = 0 ; i < nRowsToDelete ; i++) {
      rowsToDelete[i] = i + initialNumberRows_;
    }
    si_->deleteRows(nRowsToDelete, rowsToDelete);
    delete [] rowsToDelete;
    numrows_ = si_->getNumRows() ;
  }

  if (colLower_) {
    si_->setColLower(colLower_);
  }

  if (colUpper_) {
    si_->setColUpper(colUpper_);
  }

  if (cutoff_<COIN_DBL_MAX) {
    si_->setDblParam(OsiDualObjectiveLimit, cutoff_);
  }

  if (warm_) {
    if (si_->setWarmStart(warm_)==false) {
      throw CoinError("Fail restoring the warm start at the end of procedure",
          "restore","OaDecompositionBase::SaveSolverState") ;
    }
  }
  getCached();
}

void
OaDecompositionBase::passInMessageHandler(CoinMessageHandler * handler)
{
  int logLevel = handler_->logLevel();
  delete handler_;
  handler_=handler->clone();
  handler_->setLogLevel(logLevel);
}

/** Standard cut generation methods. */
void
OaDecompositionBase::generateCuts(const OsiSolverInterface &si,  OsiCuts & cs,
    const CglTreeInfo info) {
  if (nlp_ == NULL) {
    throw CoinError("Error in cut generator for outer approximation no NLP ipopt assigned", "generateCuts", "OaDecompositionBase");
  }

  // babInfo is used to communicate with the b-and-b solver (Cbc or Bcp).
  BabInfo * babInfo = dynamic_cast<BabInfo *> (si.getAuxiliaryInfo());
  assert(babInfo);
  assert(babInfo->babPtr());
  numSols_ = babInfo->babPtr()->model().getSolutionCount ();
  CglTreeInfo info_copy = info;
  const CbcNode * node = babInfo->babPtr()->model().currentNode();
  info_copy.level = (node == NULL) ? 0 : babInfo->babPtr()->model().currentNode()->depth();
  if(babInfo->hasSolution()) numSols_ ++;
  if (babInfo)
    if (!babInfo->mipFeasible())
      return;

  //Get the continuous solution
  const double *colsol = si.getColSolution();


  vector<double> savedColLower(nlp_->getNumCols());
  CoinCopyN(nlp_->getColLower(), nlp_->getNumCols(), savedColLower());
  vector<double> savedColUpper(nlp_->getNumCols());
  CoinCopyN(nlp_->getColUpper(), nlp_->getNumCols(), savedColUpper());


  OsiBranchingInformation brInfo(nlp_, false);
  brInfo.solution_ = colsol;
  //Check integer infeasibility
  bool isInteger = integerFeasible(*nlp_, brInfo, parameters_.cbcIntegerTolerance_,
                              objects_, nObjects_);


  //Check nodeNumber if it did not change scan savedCuts_ if one is violated force it and exit
  int nodeNumber = babInfo->babPtr()->model().getNodeCount();
  if(nodeNumber == currentNodeNumber_){
#ifdef OA_DEBUG
    printf("OA decomposition recalled from the same node!\n");
#endif
    int numCuts = savedCuts_.sizeRowCuts();
    for(int i = 0 ; i < numCuts ; i++){
       //Check if cuts off solution
       if(savedCuts_.rowCut(i).violated(colsol) > 0.){
#ifdef OA_DEBUG
         printf("A violated saved cut has been found\n");
#endif
         savedCuts_.rowCut(i).setEffectiveness(9.99e99);
         cs.insert(savedCuts_.rowCut(i));
         savedCuts_.eraseRowCut(i);
         return;
         i--; numCuts--;
       }
    }
  }
  else {
    currentNodeNumber_ = nodeNumber;
    savedCuts_.dumpCuts();
  } 
         
  if (!isInteger) {
    if (!doLocalSearch(babInfo))//create sub mip solver.
      return;
  }

  //get the current cutoff
  double cutoff;
  si.getDblParam(OsiDualObjectiveLimit, cutoff);

  // Save solvers state if needed

  solverManip * lpManip = NULL;
  if (lp_ != NULL) {
      assert(lp_ == &si);
      lpManip = new solverManip(lp_, true, leaveSiUnchanged_, true, true);
  }
  else {
    lpManip = new solverManip(si);
  }
  lpManip->setObjects(objects_, nObjects_);

  double milpBound = performOa(cs, *lpManip, babInfo, cutoff, info_copy);

  if(babInfo->hasSolution()){
     babInfo->babPtr()->model().setSolutionCount (numSols_ - 1);
  }

  //Transmit the bound found by the milp
  {
    if (milpBound>-1e100)
    {
      // Also store into solver
      if (babInfo)
        babInfo->setMipBound(milpBound);
    }
  }  //Clean everything :

  //  Reset the two solvers
  if (leaveSiUnchanged_)
    lpManip->restore();
  delete lpManip;

  nlp_->setColLower(savedColLower());
  nlp_->setColUpper(savedColUpper());

  return;
}

void
OaDecompositionBase::solverManip::getCached(){
  numrows_ = si_->getNumRows();
  numcols_ = si_->getNumCols();
  siColLower_ = si_->getColLower();
  siColUpper_ = si_->getColUpper();
}


/** Do update after an nlp has been solved*/
bool
OaDecompositionBase::post_nlp_solve(BabInfo * babInfo, double cutoff) const{
  nSolve_++;
  bool return_value = false;
  if (nlp_->isProvenOptimal()) {
    handler_->message(FEASIBLE_NLP, messages_)
    <<nlp_->getIterationCount()
    <<nlp_->getObjValue()<<CoinMessageEol;

#ifdef OA_DEBUG
    const double * colsol2 = nlp_->getColSolution();
    debug_.checkInteger(*nlp_,std::cerr);
#endif

    if ((nlp_->getObjValue() < cutoff) ) {
      handler_->message(UPDATE_UB, messages_)
      <<nlp_->getObjValue()
      <<CoinCpuTime()-timeBegin_
      <<CoinMessageEol;

      return_value = true;
      // Also pass it to solver
      assert(babInfo);
      if (babInfo) {
        int numcols = nlp_->getNumCols();
        double * lpSolution = new double[numcols + 1];
        CoinCopyN(nlp_->getColSolution(), numcols, lpSolution);
        lpSolution[numcols] = nlp_->getObjValue();
        babInfo->setSolution(lpSolution,
            numcols + 1, lpSolution[numcols]);
        delete [] lpSolution;
      }
    }
  }
  else if (nlp_->isAbandoned() || nlp_->isIterationLimitReached()) {
    (*handler_)<<"Unsolved NLP... exit"<<CoinMessageEol;
  }
  else {
    handler_->message(INFEASIBLE_NLP, messages_)
    <<nlp_->getIterationCount()
    <<CoinMessageEol;
  }
  return return_value;
}

void 
OaDecompositionBase::setupMipSolver(BabSetupBase &b, const std::string & prefix){


}

#ifdef OA_DEBUG
bool
OaDecompositionBase::OaDebug::checkInteger(const OsiSolverInterface &nlp, 
                                           std::ostream & os) const {
   const double * colsol = nlp.getColSolution();
   int numcols = nlp.getNumCols();
  for (int i = 0 ; i < numcols ; i++) {
    if (nlp.isInteger(i)) {
      if (fabs(colsol[i]) - floor(colsol[i] + 0.5) >
          1e-07) {
        std::cerr<<"Integer infeasible point (should not be), integer infeasibility for variable "<<i
        <<" is, "<<fabs(colsol[i] - floor(colsol[i] + 0.5))<<std::endl;
      }
    }
    return true;
  }

}

void
OaDecompositionBase::OaDebug::printEndOfProcedureDebugMessage(const OsiCuts &cs,
    bool foundSolution,
    double solValue,
    double milpBound,
    bool isInteger,
    bool feasible,
    std::ostream & os) const{
  std::cout<<"------------------------------------------------------------------"
  <<std::endl;
  std::cout<<"OA procedure finished"<<std::endl;
  std::cout<<"Generated "<<cs.sizeRowCuts()<<std::endl;
  if (foundSolution)
    std::cout <<"Found NLP-integer feasible solution of  value : "<<solValue<<std::endl;
  std::cout<<"Current MILP lower bound is : "<<milpBound<<std::endl;
  std::cout<<"-------------------------------------------------------------------"<<std::endl;
  std::cout<<"Stopped because : isInteger "<<isInteger<<", feasible "<<feasible<<std::endl<<std::endl;

}



#endif
}/* End namespace Bonmin. */

