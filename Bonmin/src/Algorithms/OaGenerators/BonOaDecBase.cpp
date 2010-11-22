// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
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
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h" 
#define CHECK_CPX_STAT(a,b) if(b) throw CoinError("Error in CPLEX call",__FILE__,a);

#endif
#include "OsiAuxInfo.hpp"

//The following two are to interupt the solution of sub-mip through CTRL-C
extern CbcModel * OAModel;

namespace Bonmin
{

  OaDecompositionBase::OaDecompositionBase
  (OsiTMINLPInterface * nlp,
   OsiSolverInterface * si,
   CbcStrategy * strategy,
   double cbcCutoffIncrement,
   double cbcIntegerTolerance,
   bool leaveSiUnchanged
  )
      :
      CglCutGenerator(),
      nlp_(nlp),
      nSolve_(0),
      lp_(si),
      objects_(NULL),
      nObjects_(0),
      nLocalSearch_(0),
      handler_(NULL),
      leaveSiUnchanged_(leaveSiUnchanged),
      reassignLpsolver_(false),
      timeBegin_(0),
      parameters_()
  {
    handler_ = new CoinMessageHandler();
    handler_ -> setLogLevel(2);
    messages_ = OaMessages();
    if (strategy)
      parameters_.setStrategy(*strategy);
    timeBegin_ = CoinCpuTime();
    parameters_.cbcCutoffIncrement_  = cbcCutoffIncrement;
    parameters_.cbcIntegerTolerance_ = cbcIntegerTolerance;
  }
  OaDecompositionBase::OaDecompositionBase(BabSetupBase &b, bool leaveSiUnchanged,
      bool reassignLpsolver):
      CglCutGenerator(),
      nlp_(b.nonlinearSolver()),
      lp_(b.continuousSolver()),
      objects_(NULL),
      nObjects_(0),
      nLocalSearch_(0),
      handler_(NULL),
      leaveSiUnchanged_(leaveSiUnchanged),
      reassignLpsolver_(reassignLpsolver),
      timeBegin_(0),
      parameters_()
  {
    handler_ = new CoinMessageHandler();
    int logLevel;
    b.options()->GetIntegerValue("oa_log_level",logLevel,"bonmin.");
    b.options()->GetNumericValue("oa_log_frequency",parameters_.logFrequency_,"bonmin.");

    handler_ -> setLogLevel(logLevel);

    messages_ = OaMessages();
    timeBegin_ = CoinCpuTime();
    b.options()->GetIntegerValue("milp_log_level",parameters_.subMilpLogLevel_,"bonmin.");
    b.options()->GetNumericValue("cutoff_decr",parameters_.cbcCutoffIncrement_,"bonmin.");
    b.options()->GetNumericValue("integer_tolerance",parameters_.cbcIntegerTolerance_,"bonmin.");
    int ivalue;
    b.options()->GetEnumValue("add_only_violated_oa", ivalue,"bonmin.");
    parameters_.addOnlyViolated_ = ivalue;
    b.options()->GetEnumValue("oa_cuts_scope", ivalue,"bonmin.");
    parameters_.global_ = ivalue;
  }

  OaDecompositionBase::OaDecompositionBase
  (const OaDecompositionBase & other)
      :
      CglCutGenerator(other),
      nlp_(other.nlp_),
      lp_(other.lp_),
      objects_(other.objects_),
      nObjects_(other.nObjects_),
      nLocalSearch_(0),
      messages_(other.messages_),
      leaveSiUnchanged_(other.leaveSiUnchanged_),
      reassignLpsolver_(other.reassignLpsolver_),
      timeBegin_(other.timeBegin_),
      parameters_(other.parameters_)
  {
    //timeBegin_ = CoinCpuTime();
    handler_ = other.handler_->clone();
  }
/// Constructor with default values for parameters
  OaDecompositionBase::Parameters::Parameters():
      global_(true),
      addOnlyViolated_(false),
      cbcCutoffIncrement_(1e-06),
      cbcIntegerTolerance_(1e-05),
      localSearchNodeLimit_(0),
      maxLocalSearchPerNode_(0),
      maxLocalSearch_(0),
      maxLocalSearchTime_(3600),
      subMilpLogLevel_(0),
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
      localSearchNodeLimit_(other.localSearchNodeLimit_),
      maxLocalSearchPerNode_(other.maxLocalSearchPerNode_),
      maxLocalSearch_(other.maxLocalSearch_),
      maxLocalSearchTime_(other.maxLocalSearchTime_),
      subMilpLogLevel_(other.subMilpLogLevel_),
      logFrequency_(other.logFrequency_),
      strategy_(NULL)
  {
    if (other.strategy_)
      strategy_ = other.strategy_->clone();
  }


  /** Constructor */
  OaDecompositionBase::SubMipSolver::SubMipSolver(OsiSolverInterface * lp,
      const CbcStrategy * strategy):
      lp_(lp),
      clp_(NULL),
      cpx_(NULL),
      cbc_(NULL),
      lowBound_(-COIN_DBL_MAX),
      optimal_(false),
      integerSolution_(NULL),
      strategy_(NULL)
  {
    clp_ = (lp_ == NULL)? NULL :
        dynamic_cast<OsiClpSolverInterface *>(lp_);
#ifdef COIN_HAS_CPX
    cpx_ = (lp_ == NULL)? NULL :
        dynamic_cast<OsiCpxSolverInterface *>(lp_);
#endif
    if (strategy) strategy_ = strategy->clone();
  }
  OaDecompositionBase::SubMipSolver::~SubMipSolver()
  {
    if (strategy_) delete strategy_;
    if (integerSolution_) delete [] integerSolution_;
    if (cbc_) delete cbc_;
  }

  /** Assign lp solver. */
  void
  OaDecompositionBase::SubMipSolver::setLpSolver(OsiSolverInterface * lp)
  {
    lp_ = lp;
    clp_ = (lp_ == NULL) ? NULL :
        dynamic_cast<OsiClpSolverInterface *>(lp_);
#ifdef COIN_HAS_CPX
    cpx_ = (lp_ == NULL) ? NULL :
        dynamic_cast<OsiCpxSolverInterface *>(lp_);
#endif
    lowBound_ = -COIN_DBL_MAX;
    optimal_ = false;
    if (integerSolution_) {
      delete [] integerSolution_;
      integerSolution_ = NULL;
    }
  }



  void
  OaDecompositionBase::SubMipSolver::performLocalSearch(double cutoff, int loglevel, double maxTime,
      int maxNodes)
  {
    if (clp_) {
      if (!strategy_)
        strategy_ = new CbcStrategyDefault(1,0,0, loglevel);

      OsiBabSolver empty;
      if (cbc_) delete cbc_;
      OAModel = cbc_ = new CbcModel(*clp_);
      cbc_->solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(cbc_->messagesPointer()->source_,"OaCbc");

      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      clp_->resolve();
      cbc_->setStrategy(*strategy_);
      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      cbc_->setMaximumNodes(maxNodes);
      cbc_->setMaximumSeconds(maxTime);
      cbc_->setCutoff(cutoff);
      cbc_->branchAndBound();
      OAModel = NULL;
      lowBound_ = cbc_->getBestPossibleObjValue();

      if (cbc_->isProvenOptimal() || cbc_->isProvenInfeasible())
        optimal_ = true;
      else optimal_ = false;

      if (cbc_->getSolutionCount()) {
        if (!integerSolution_)
          integerSolution_ = new double[lp_->getNumCols()];
        CoinCopyN(cbc_->bestSolution(), lp_->getNumCols(), integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
      nodeCount_ = cbc_->getNodeCount();
      iterationCount_ = cbc_->getIterationCount();
    }
    else {
      lp_->messageHandler()->setLogLevel(loglevel);
#ifdef COIN_HAS_CPX
      if (cpx_) {
        CPXENVptr env = cpx_->getEnvironmentPtr();
        CPXsetintparam(env, CPX_PARAM_NODELIM, maxNodes);
        CPXsetdblparam(env, CPX_PARAM_TILIM, maxTime);
        CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);
        //CpxModel = cpx_;
      }
      else
#endif 
     {
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
    }

#ifdef COIN_HAS_CPX
    if (cpx_) {
      cpx_->switchToMIP();
      //CpxModel = NULL;
      CPXENVptr env = cpx_->getEnvironmentPtr();
      CPXLPptr cpxlp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

      int status = CPXmipopt(env,cpxlp);
      CHECK_CPX_STAT("mipopt",status)

      int stat = CPXgetstat( env, cpxlp);
//#define FOR_HASSAN
#ifdef FOR_HASSAN
      std::cout << "STATUS = "<<stat<<"\n";
      if(stat==119){
          cpx_->writeMpsNative("OA.mps", NULL, NULL, 1);
          throw CoinError("Ok program stoped by me !","OaDecompositionBase::SubMipSolver","performLocalSearch");
      }
#endif
      bool infeasible = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_ABORT_INFEAS) || (stat == CPXMIP_TIME_LIM_INFEAS) || (stat == CPXMIP_NODE_LIM_INFEAS) || (stat == CPXMIP_FAIL_INFEAS)
                        || (stat == CPXMIP_MEM_LIM_INFEAS) || (stat == CPXMIP_INForUNBD);
      optimal_ = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL); 
      nodeCount_ = CPXgetnodecnt(env , cpxlp);
      iterationCount_ = CPXgetmipitcnt(env , cpxlp);
      status = CPXgetbestobjval(env, cpxlp, &lowBound_);
      CHECK_CPX_STAT("getbestobjval",status)
       
      if(!infeasible){
         if(!integerSolution_){
           integerSolution_ = new double[lp_->getNumCols()];
         }
         CPXgetmipx(env, cpxlp, integerSolution_, 0, lp_->getNumCols() -1);
         CHECK_CPX_STAT("getmipx",status)
      }
      else {
        if (integerSolution_) {
          delete [] integerSolution_;
          integerSolution_ = NULL;
        }
      }
    }
    else {
#endif
    lp_->branchAndBound();

   optimal_ = lp_->isProvenOptimal();

    if (lp_->getFractionalIndices().size() == 0) {
      if (!integerSolution_)
        integerSolution_ = new double[lp_->getNumCols()];
      CoinCopyN(lp_->getColSolution(), lp_->getNumCols() , integerSolution_);
    }
    else if (integerSolution_) {
      delete [] integerSolution_;
      integerSolution_ = NULL;
    }
  }
#ifdef COIN_HAS_CPX
  }
#endif
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
    cutoff_(COIN_DBL_MAX),
    deleteSolver_(false),
    objects_(NULL),
    nObjects_(0),
    integerTolerance_(1e-08)
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
    cutoff_(COIN_DBL_MAX),
    deleteSolver_(true),
    objects_(NULL),
    nObjects_(0),
    integerTolerance_(1e-08)
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
    int nRowsToDelete = numrows_ - initialNumberRows_;
    int * rowsToDelete = new int[nRowsToDelete];
    for (int i = 0 ; i < nRowsToDelete ; i++) {
      rowsToDelete[i] = i + initialNumberRows_;
    }
    si_->deleteRows(nRowsToDelete, rowsToDelete);
    delete [] rowsToDelete;
    numrows_ -= nRowsToDelete;
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

/// Check for integer feasibility of a solution return 1 if it is
bool
OaDecompositionBase::solverManip::integerFeasible(const OsiBranchingInformation& info) const
{
  if (objects_) {
    int dummy;
    for (int i = 0 ; i < nObjects_ ; i++) {
      if (objects_[i]->infeasibility(&info, dummy) > 0.0) return false;
    }
  }
  else {
    const double * sol = info.solution_;
    int numcols = si_->getNumCols();
    for (int i = 0 ; i < numcols ; i++) {
      if (si_->isInteger(i)) {
        if (fabs(sol[i] - floor(sol[i] + 0.5)) >
            integerTolerance_) {
          return false;
        }
      }
    }
  }
  return true;
}

/** Fix integer variables in si to their values in colsol.
*/
void
OaDecompositionBase::solverManip::fixIntegers(const OsiBranchingInformation& info)
{
  if (objects_) {
    for (int i = 0 ; i < nObjects_ ; i++) {
      if (objects_[i]->feasibleRegion(si_, &info));
    }
  }
  else {
    const double * colsol = info.solution_;
    for (int i = 0; i < numcols_; i++) {
      if (si_->isInteger(i)) {
        double  value =  colsol[i];
        if (fabs(value - floor(value+0.5)) > integerTolerance_) {
          std::stringstream stream;
          stream<<"Error not integer valued solution"<<std::endl;
          stream<<"---------------- x["<<i<<"] = "<<value<<std::endl;
          throw CoinError(stream.str(),"fixIntegers","OaDecompositionBase::solverManip");
        }
        value = floor(value+0.5);
        value = std::max(colLower_[i],value);
        value = std::min(value, colUpper_[i]);

        if (fabs(value) > 1e10) {
          std::stringstream stream;
          stream<<"Can not fix variable in nlp because it has too big a value ("<<value
          <<") at optimium of LP relaxation. You should try running the problem with B-BB"<<std::endl;
          throw CoinError(stream.str(),
              "fixIntegers","OaDecompositionBase::solverManip") ;
        }
#ifdef OA_DEBUG
        //         printf("xx %d at %g (bounds %g, %g)",i,value,nlp_->getColLower()[i],
        //                nlp_->getColUpper()[i]);
        std::cout<<(int)value;
#endif
        si_->setColLower(i,value);
        si_->setColUpper(i,value);
      }
    }
#ifdef OA_DEBUG
    std::cout<<std::endl;
#endif
  }
}
/** Check if solution in solver is the same as colsol on integer variables. */
bool
OaDecompositionBase::solverManip::isDifferentOnIntegers(const double * colsol)
{
  return isDifferentOnIntegers(colsol, si_->getColSolution());
}

/** Check if two solutions are the same on integer variables. */
bool
OaDecompositionBase::solverManip::isDifferentOnIntegers(const double * colsol, const double *otherSol)
{
  if (objects_) {
    for (int i = 0 ; i < nObjects_ ; i++) {
      OsiObject * obj = objects_[i];
      int colnum = obj->columnNumber();
      if (colnum >= 0) {//Variable branching object
        if (fabs(otherSol[colnum] - colsol[colnum]) > 100*integerTolerance_) {
          return true;
        }
      }
      else {//It is a sos branching object
        OsiSOS * sos = dynamic_cast<OsiSOS *>(obj);
        assert(sos);
        const int * members = sos->members();
        int end = sos->numberMembers();
        for (int k = 0 ; k < end ; k++) {
          if (fabs(otherSol[members[k]] - colsol[members[k]]) > 0.001) {
            return true;
          }
        }
      }
    }
  }
  else {
    for (int i = 0; i < numcols_ ; i++) {
      if (si_->isInteger(i) && fabs(otherSol[i] - colsol[i])>0.001)
        return true;
    }
  }
  return false;
}

/** Clone the state of another solver (bounds, cutoff, basis).*/
void
OaDecompositionBase::solverManip::cloneOther(const OsiSolverInterface &si)
{
  //Install current active cuts into local solver
  int numberCutsToAdd = si.getNumRows();
  numberCutsToAdd -= numrows_;
  if (numberCutsToAdd > 0)//Have to install some cuts
  {
    CoinPackedVector * * addCuts = new CoinPackedVector *[numberCutsToAdd];
    for (int i = 0 ; i < numberCutsToAdd ; i++)
    {
      addCuts[i] = new CoinPackedVector;
    }
    //Get the current matrix and fill the addCuts
    const CoinPackedMatrix * mat = si.getMatrixByCol();
    const CoinBigIndex * start = mat->getVectorStarts();
    const int * length = mat->getVectorLengths();
    const double * elements = mat->getElements();
    const int * indices = mat->getIndices();
    for (int i = 0 ; i < numcols_ ; i++)
      for (int k = start[i] ; k < start[i] + length[i] ; k++)
      {
        if (indices[k] >= numrows_) {
          addCuts[ indices[k] - numrows_ ]->insert(i, elements[k]);
        }
      }
    si_->addRows(numberCutsToAdd, (const CoinPackedVectorBase * const *) addCuts, si.getRowLower() + numrows_,
        si.getRowUpper() + numrows_);
    for (int i = 0 ; i < numberCutsToAdd ; i++){
      delete addCuts[i];
    }
    delete [] addCuts;
  }
  else if (numberCutsToAdd < 0) { //Oups some error
    throw CoinError("Internal error number of cuts wrong",
        "generateCuts","OACutGenerator2");
  }

  si_->setColLower(si.getColLower());
  si_->setColUpper(si.getColUpper());
  //Install basis in problem
  CoinWarmStart * warm = si.getWarmStart();
  if (si_->setWarmStart(warm)==false) {
    delete warm;
    throw CoinError("Fail installing the warm start in the subproblem",
        "generateCuts","OACutGenerator2") ;
  }
  delete warm;
  //put the cutoff
  double cutoff;
  si.getDblParam(OsiDualObjectiveLimit, cutoff);
  si_->setDblParam(OsiDualObjectiveLimit, cutoff);
  si_->messageHandler()->setLogLevel(0);
  si_->resolve();

  numrows_ = si.getNumRows();
#ifdef OA_DEBUG

  std::cout<<"Resolve with hotstart :"<<std::endl
  <<"number of iterations(should be 0) : "<<lp_->getIterationCount()<<std::endl
  <<"Objective value and diff to original : "<<lp_->getObjValue()<<", "
  <<fabs(si_->getObjValue() - si.getObjValue())<<std::endl;
  for (int i = 0 ; i <= numcols ; i++) {
    if (fabs(si.getColSolution()[i]-si_->getColSolution()[i])>1e-08) {
      std::cout<<"Diff between solution at node and solution with local solver : "<<fabs(si.getColSolution()[i]-si_->getColSolution()[i])<<std::endl;
    }
  }
#endif

}


/** Install cuts in solver. */
void
OaDecompositionBase::solverManip::installCuts(const OsiCuts& cs, int numberCuts)
{
  int numberCutsBefore = cs.sizeRowCuts() - numberCuts;

  CoinWarmStartBasis * basis
  = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart()) ;
  assert(basis != NULL); // make sure not volume
  basis->resize(numrows_ + numberCuts,numcols_) ;
  for (int i = 0 ; i < numberCuts ; i++) {
    basis->setArtifStatus(numrows_ + i,
        CoinWarmStartBasis::basic) ;
  }

  const OsiRowCut ** addCuts = new const OsiRowCut * [numberCuts] ;
  for (int i = 0 ; i < numberCuts ; i++) {
    addCuts[i] = &cs.rowCut(i + numberCutsBefore) ;
  }
  si_->applyRowCuts(numberCuts,addCuts) ;
  numrows_ += numberCuts;
  delete [] addCuts ;
  if (si_->setWarmStart(basis) == false) {
    delete basis;
    throw CoinError("Fail setWarmStart() after cut installation.",
        "generateCuts","OACutGenerator2") ;
  }
  delete basis;
}

/** Standard cut generation methods. */
void
OaDecompositionBase::generateCuts(const OsiSolverInterface &si,  OsiCuts & cs,
    const CglTreeInfo info) const
{
  if (nlp_ == NULL) {
    throw CoinError("Error in cut generator for outer approximation no NLP ipopt assigned", "generateCuts", "OaDecompositionBase");
  }

  // babInfo is used to communicate with the b-and-b solver (Cbc or Bcp).
  OsiBabSolver * babInfo = dynamic_cast<OsiBabSolver *> (si.getAuxiliaryInfo());
  if (babInfo)
    if (!babInfo->mipFeasible())
      return;

  //Get the continuous solution
  const double *colsol = si.getColSolution();


  solverManip nlpManip(nlp_, false, false, true, false, false);
  nlpManip.setIntegerTolerance(parameters_.cbcIntegerTolerance_);
  nlpManip.setObjects(objects_, nObjects_);
  OsiBranchingInformation brInfo(nlp_, false);
  brInfo.solution_ = colsol;
  //Check integer infeasibility
  bool isInteger = nlpManip.integerFeasible(brInfo);

  SubMipSolver * subMip = NULL;

  if (!isInteger) {
    if (doLocalSearch())//create sub mip solver.
    {
      subMip = new SubMipSolver(lp_, parameters_.strategy());
    }
    else {
      return;
    }
  }


  //If we are going to modify things copy current information to restore it in the end


  //get the current cutoff
  double cutoff;
  si.getDblParam(OsiDualObjectiveLimit, cutoff);

  // Save solvers state if needed

  solverManip * lpManip = NULL;
  if (lp_ != NULL) {
    if (lp_!=&si) {
      lpManip = new solverManip(lp_, true, false, false, true, true);
      lpManip->cloneOther(si);
    }
    else {
#if 0
      throw CoinError("Not allowed to modify si in a cutGenerator",
          "OACutGenerator2","generateCuts");
#else
      lpManip = new solverManip(lp_, true, leaveSiUnchanged_, true, true);
#endif
    }
  }
  else {
    lpManip = new solverManip(si);
  }
  lpManip->setObjects(objects_, nObjects_);
  lpManip->setIntegerTolerance(parameters_.cbcIntegerTolerance_);
  double milpBound = performOa(cs, nlpManip, *lpManip, subMip, babInfo, cutoff);

  //Transmit the bound found by the milp
  {
    if (milpBound>-1e100)
    {
      // Also store into solver
      if (babInfo)
        babInfo->setMipBound(milpBound);
    }
  }  //Clean everything :

  //free subMip
  if (subMip!= NULL) {
    delete subMip;
    subMip = NULL;
  }

  //  Reset the two solvers
  if (leaveSiUnchanged_)
    lpManip->restore();
  delete lpManip;
  nlpManip.restore();
  return;
}

void
OaDecompositionBase::solverManip::getCached()
{
  numrows_ = si_->getNumRows();
  numcols_ = si_->getNumCols();
  siColLower_ = si_->getColLower();
  siColUpper_ = si_->getColUpper();
}


/** Solve the nlp and do output return true if feasible*/
bool
OaDecompositionBase::solveNlp(OsiBabSolver * babInfo, double cutoff) const
{
  nSolve_++;
  nlp_->resolve();
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
      if (babInfo) {
        int numcols = nlp_->getNumCols();
        double * lpSolution = new double[numcols + 1];
        CoinCopyN(nlp_->getColSolution(), numcols, lpSolution);
        lpSolution[numcols] = nlp_->getObjValue();
        babInfo->setSolution(lpSolution,
            numcols + 1, lpSolution[numcols]);
        delete [] lpSolution;
      }
      else {
        printf("No auxiliary info in nlp solve!\n");
        throw -1;
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



#ifdef OA_DEBUG
bool
OaDecompositionBase::OaDebug::checkInteger(const OsiSolverInterface &nlp, std::ostream & os) const
{
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
    std::ostream & os) const
{
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

