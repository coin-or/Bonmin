// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 05/26/2005
//#define OA_DEBUG

#include "BonminConfig.h"

#include "BonOACutGenerator2.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#include "OsiAuxInfo.hpp"


extern CbcModel * OAModel;
namespace Bonmin{
/// Default constructor
OACutGenerator2::OACutGenerator2():
    CglCutGenerator(),
    nlp_(NULL),
    nSolve_(0),
    si_(NULL),
    cbcCutoffIncrement_(1e-06),
    cbcIntegerTolerance_(1e-05),
    localSearchNodeLimit_(0),
    maxLocalSearchPerNode_(0),
    maxLocalSearch_(0),
    maxLocalSearchTime_(3600),
    nLocalSearch_(0),
    solveAuxiliaryProblem_(0),
    handler_(NULL),
    subMilpLogLevel_(0),
    leaveSiUnchanged_(0),
    strategy_(NULL),
    timeBegin_(0.),
    logFrequency_(1000.)
{
  handler_ = new CoinMessageHandler();
  handler_ -> setLogLevel(2);
  messages_ = OaMessages();
  timeBegin_ = CoinCpuTime();
}



OACutGenerator2::OACutGenerator2
(OsiTMINLPInterface * nlp,
 OsiSolverInterface * si,
 CbcStrategy * strategy,
 double cbcCutoffIncrement,
 double cbcIntegerTolerance,
 bool solveAuxiliaryProblem,
 bool leaveSiUnchanged
)
    :
    CglCutGenerator(),
    nlp_(nlp),
    nSolve_(0),
    si_(si),
    cbcCutoffIncrement_(cbcCutoffIncrement),
    cbcIntegerTolerance_(cbcIntegerTolerance),
    localSearchNodeLimit_(0),
    maxLocalSearchPerNode_(0),
    maxLocalSearch_(0),
    maxLocalSearchTime_(3600),
    nLocalSearch_(0),
    solveAuxiliaryProblem_(solveAuxiliaryProblem),
    handler_(NULL),
    subMilpLogLevel_(0),
    leaveSiUnchanged_(leaveSiUnchanged),
    strategy_(NULL),
    timeBegin_(0),
    logFrequency_(1000.)
{
  handler_ = new CoinMessageHandler();
  handler_ -> setLogLevel(2);
  messages_ = OaMessages();
  if(strategy)
    strategy_ = strategy->clone();
  timeBegin_ = CoinCpuTime();
}

OACutGenerator2::~OACutGenerator2()
{
  delete handler_;
  if(strategy_)
    delete strategy_;
}
/// Assign an OsiTMINLPInterface
void
OACutGenerator2::assignNlpInterface(OsiTMINLPInterface * nlp)
{
  nlp_ = nlp;
}

/// Assign an OsiTMINLPInterface
void
OACutGenerator2::assignLpInterface(OsiSolverInterface * si)
{
  si_ = si;
  if(maxLocalSearch_>0) {
    setTheNodeLimit();
  }
}

double OACutGenerator2::siBestObj(CbcModel * model) const
{
  if(model == NULL) {
    //Check solver name to see if local searches can be performed
    std::string solverName;
    si_->getStrParam(OsiSolverName,solverName);
    //for cbc needs to get the first three letter of solver name
    std::string shortSolverName(solverName,0,3);
    if(shortSolverName == "cbc") {
      throw CoinError("OsiCbc is not supported for doing local searches use OsiClpSolverInterface instead",
          "OACutGenerator2","setTheNodeLimit");
    }
    else if(solverName == "cplex") {
#ifdef COIN_HAS_CPX
      OsiCpxSolverInterface * cpx = dynamic_cast<OsiCpxSolverInterface *>(si_);
      double value;
      int status = CPXgetbestobjval(cpx->getEnvironmentPtr(),cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &value);
      if(status)
        throw CoinError("Error in getting CPLEX best bound","OACutGenerator2","siBestObj");
      return value;
#else

      throw CoinError("You have to define COIN_HAS_CPX at compilation to be able to use cplex for local searches",
          "OACutGenerator2","siBestObj");
#endif

    }
    else {
      std::string mesg="Unsuported solver";
      mesg += solverName;
      mesg +=" for local searches, you should use Cbc or Cplex";
      throw CoinError(mesg,
          "OACutGenerator2","assignLpInterface");
    }
  }
  else
    return model->getBestPossibleObjValue();
}
void OACutGenerator2::setTheNodeLimit()
{
  //Check solver name to see if local searches can be performed
  std::string solverName;
  si_->getStrParam(OsiSolverName,solverName);

  //for cbc needs to get the first three letter of solver name
  std::string shortSolverName(solverName,0,3);
  if(shortSolverName == "cbc") {
    throw CoinError("OsiCbc is not supported for doing local searches use OsiClpSolverInterface instead",
        "OACutGenerator2","setTheNodeLimit");
  }
  else if(solverName == "cplex") {
#ifdef COIN_HAS_CPX
    OsiCpxSolverInterface * cpx = dynamic_cast<OsiCpxSolverInterface *>(si_);
    CPXsetintparam(cpx->getEnvironmentPtr(),CPX_PARAM_NODELIM, localSearchNodeLimit_);
#else

    throw CoinError("You have to define COIN_HAS_CPX at compilation to be able to use cplex for local searches",
        "OACutGenerator2","setTheNodeLimit");
#endif

  }
  else if(solverName == "clp") {
    //do nothing will set up node limit when creating a CbcModel
  }
  else {
    std::string mesg="Unsuported solver";
    mesg += solverName;
    mesg +=" for local searches, you should use Cbc or Cplex";
    throw CoinError(mesg,
        "OACutGenerator2","setTheNodeLimit");
  }
}


void OACutGenerator2::setTimeLimit(double time) const
{
  //Check solver name to see if local searches can be performed
  std::string solverName;
  si_->getStrParam(OsiSolverName,solverName);
  //for cbc needs to get the first three letter of solver name
  std::string shortSolverName(solverName,0,3);
  if(shortSolverName == "cbc") {
    throw CoinError("OsiCbc is not supported for doing local searches use OsiClpSolverInterface instead",
        "OACutGenerator2","setTheNodeLimit");
  }
  else if(solverName == "cplex") {
#ifdef COIN_HAS_CPX
    OsiCpxSolverInterface * cpx = dynamic_cast<OsiCpxSolverInterface *>(si_);
    CPXsetdblparam(cpx->getEnvironmentPtr(),CPX_PARAM_TILIM, time);
#else

    throw CoinError("You have to define COIN_HAS_CPX at compilation to be able to use cplex for local searches",
        "OACutGenerator2","setTheNodeLimit");
#endif

  }
  else if(solverName == "clp") {
    //do nothing will set up node limit when creating a CbcModel
  }
  else {
    std::string mesg="Unsuported solver";
    mesg += solverName;
    mesg +=" for local searches, you should use Cbc or Cplex";
    throw CoinError(mesg,
        "OACutGenerator2","setTheNodeLimit");
  }
}

void OACutGenerator2::setCutoff(double bestKnown) const
{
  //Check solver name to see if local searches can be performed
  std::string solverName;
  si_->getStrParam(OsiSolverName,solverName);

  //for cbc needs to get the first three letter of solver name
  std::string shortSolverName(solverName,0,3);
  if(shortSolverName == "cbc") {
    throw CoinError("OsiCbc is not supported for doing local searches use OsiClpSolverInterface instead",
        "OACutGenerator2","setTheNodeLimit");
  }
  else if(solverName == "cplex") {
#ifdef COIN_HAS_CPX
    OsiCpxSolverInterface * cpx = dynamic_cast<OsiCpxSolverInterface *>(si_);
    CPXsetdblparam(cpx->getEnvironmentPtr(),CPX_PARAM_CUTUP, bestKnown);
#else

    throw CoinError("You have to define COIN_HAS_CPX at compilation to be able to use cplex for local searches",
        "OACutGenerator2","setTheNodeLimit");
#endif

  }
  else if(solverName == "clp") {
    //do nothing will set up node limit when creating a CbcModel
  }
  else {
    std::string mesg="Unsuported solver";
    mesg += solverName;
    mesg +=" for local searches, you should use Cbc or Cplex";
    throw CoinError(mesg,
        "OACutGenerator2","setTheNodeLimit");
  }
}/// cut generation method
void
OACutGenerator2::generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
    const CglTreeInfo info) const
{
  CbcStrategy * strategy = strategy_;
  double lastPeriodicLog= CoinCpuTime();

  if(nlp_ == NULL) {
    std::cerr<<"Error in cut generator for outer approximation no NLP ipopt assigned"<<std::endl;
    throw -1;
  }

  //This node may already be fathomed by a previous call to this
  // this information is stored in babInfo
  OsiBabSolver * babInfo = dynamic_cast<OsiBabSolver *> (si.getAuxiliaryInfo());
#if 1
  //if algorithms would converge this should never happen and it seems very dangerous if somebody forgot to reset the bound
  if(babInfo)
    if(!babInfo->mipFeasible())
      return;

#endif

  const int numcols = nlp_->getNumCols();

  //Get the continuous solution
  const double *colsol = si.getColSolution();

  bool doLocalSearch = 0;
  bool isInteger = 1;
  //   nlp_->setMilpBound(-1e100);
  //   nlp_->setMilpFeasible(true);


  //Check integer infeasibility
  for(int i = 0 ; i < numcols ; i++) {
    if(nlp_->isInteger(i)) {
      if(fabs(colsol[i] - floor(colsol[i] + 0.5) ) >
          cbcIntegerTolerance_) {
        isInteger = 0;
        break;
      }
    }
  }

  if(!isInteger) {
    if(nLocalSearch_<maxLocalSearch_ &&
        localSearchNodeLimit_ > 0 &&
        CoinCpuTime() - timeBegin_ < maxLocalSearchTime_)//do a local search
    {
      doLocalSearch = 1;
    }
    else {
      // 	  nlp_->setMilpBound(-1e100);
      // 	  nlp_->setMilpFeasible(true);
      if(strategy && ! strategy_)
        delete strategy;
      return;
    }
  }


  //We are going to generate some cuts copy LP information


  //get the current cutoff
  double cutoff;
  si.getDblParam(OsiDualObjectiveLimit, cutoff);
  double saveCutoff=cutoff;

  // Save column bounds and basis
  double * saveColLb = new double[numcols];
  double * saveColUb = new double[numcols];
  CoinCopyN(nlp_->getColLower(), numcols , saveColLb);
  CoinCopyN(nlp_->getColUpper(), numcols , saveColUb);
  int originalRowNumber = (si_!=NULL) ?si_->getNumRows() : -1;
  CoinWarmStart * saveWarmStart = NULL;
  bool milpOptimal = 1;

#ifdef NO_NULL_SI
  if(si_ == NULL) {
    std::cerr<<"Error in cut generator for outer approximation no lp solver interface assigned"<<std::endl;
    throw -1;
  }
#else

  bool deleteSi = false;
  if(si_ == NULL) {
    si_ = si.clone();
    deleteSi = true;
  }
  else
#endif
  if(si_!=&si)//We may have the same lp as the one passed by the model or a local copy
  {

    //Install current active cuts into local solver
    int numberCutsToAdd = si.getNumRows();
    numberCutsToAdd -= si_->getNumRows();
    if(numberCutsToAdd > 0)//Have to install some cuts
    {
      CoinPackedVector * * addCuts = new CoinPackedVector *[numberCutsToAdd];
      for(int i = 0 ; i < numberCutsToAdd ; i++)
      {
        addCuts[i] = new CoinPackedVector;
      }
      //Get the current matrix and fill the addCuts
      const CoinPackedMatrix * mat = si.getMatrixByCol();
      const CoinBigIndex * start = mat->getVectorStarts();
      const int * length = mat->getVectorLengths();
      const double * elements = mat->getElements();
      const int * indices = mat->getIndices();
      for(int i = 0 ; i <= numcols ; i++)
        for(int k = start[i] ; k < start[i] + length[i] ; k++)
        {
          if(indices[k] >= si_->getNumRows()) {
            addCuts[ indices[k] - si_->getNumRows() ]->insert(i, elements[k]);
          }
        }
      si_->addRows(numberCutsToAdd, (const CoinPackedVectorBase * const *) addCuts, &si.getRowLower()[si_->getNumRows()],
          &si.getRowUpper()[si_->getNumRows()]);
    }
    else if (numberCutsToAdd < 0)//Oups some error
    {
      std::cerr<<"Internal error in OACutGenerator2 : number of cuts wrong"<<std::endl;
    }

    //Set the bounds on columns
    for(int i = 0 ; i <= numcols ; i++)
    {
      si_->setColBounds(i, si.getColLower()[i], si.getColUpper()[i]);
    }
    //Install basis in problem
    CoinWarmStart * warm = si.getWarmStart();
    if(si_->setWarmStart(warm)==false)
    {
      delete warm;
      throw CoinError("Fail installing the warm start in the subproblem",
          "generateCuts","OACutGenerator2") ;
    }
    delete warm;
    //put the cutoff
    si_->setDblParam(OsiDualObjectiveLimit, cutoff);
    si_->resolve();

#ifdef OA_DEBUG

    std::cout<<"Resolve with hotstart :"<<std::endl
    <<"number of iterations(should be 0) : "<<si_->getIterationCount()<<std::endl
    <<"Objective value and diff to original : "<<si_->getObjValue()<<", "
    <<fabs(si_->getObjValue() - si.getObjValue())<<std::endl;
    for(int i = 0 ; i <= numcols ; i++)
    {
      if(fabs(si.getColSolution()[i]-si_->getColSolution()[i])>1e-08) {
        std::cout<<"Diff between solution at node and solution with local solver : "<<fabs(si.getColSolution()[i]-si_->getColSolution()[i])<<std::endl;
      }
    }
#endif
  }
  else {
#ifdef NO_NULL_SI
    throw CoinError("Not allowed to modify si in a cutGenerator",
        "OACutGenerator2","generateCuts");
#else //Seems that nobody wants to allow me to do this
    if(leaveSiUnchanged_)
      saveWarmStart = si.getWarmStart();
#endif
  }

  OsiClpSolverInterface * clp = dynamic_cast<OsiClpSolverInterface *>(si_);
  CbcModel * model = NULL;//We will need the model through all the function
  double milpBound = -DBL_MAX;
  bool milpFeasible = 1;
  //double milpBound=si_->getObjValue();
  bool feasible = 1;

  if(feasible && doLocalSearch)//Perform a local search
  {
    if(clp)//Clp does not do branch-and-bound have to create a CbcModel
    {
      OsiBabSolver empty;
      model = new CbcModel(*clp); // which clones
      OAModel = model;
      model->solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(model->messagesPointer()->source_,"OaCbc");

      if(!strategy)
        strategy = new CbcStrategyDefault(1,0,0,subMilpLogLevel_);

      clp->resolve();
      model->setLogLevel(subMilpLogLevel_);
      model->solver()->messageHandler()->setLogLevel(0);
      model->setStrategy(*strategy);
      model->setMaximumNodes(localSearchNodeLimit_);
      model->setMaximumSeconds(maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());
      model->setCutoff(cutoff);
      model->branchAndBound();
      milpBound = siBestObj(model);
      if(model->isProvenOptimal())
      {
        milpOptimal = 1;
      }
      else
        milpOptimal = 0;
      feasible = milpBound < cutoff;
      milpFeasible = feasible;
      isInteger = model->getSolutionCount();
      nLocalSearch_++;

      if(milpOptimal)
        handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<model->getNodeCount()<<model->getIterationCount()<<CoinMessageEol;
      else
      {
        handler_->message(LOCAL_SEARCH_ABORT, messages_)<<model->getNodeCount()<<model->getIterationCount()<<CoinMessageEol;
        if(isInteger)
          std::cout<<"Integer feasible solution found"<<std::endl;
      }
    }
    else//use OsiSolverInterface::branchAndBound
    {
      setTimeLimit(maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());
      si_->messageHandler()->setLogLevel(subMilpLogLevel_);
      si_->branchAndBound();
      nLocalSearch_++;
      //Did we find any integer feasible solution
      isInteger = si_->getFractionalIndices().size() == 0;
      milpBound = siBestObj();
      feasible = milpBound < cutoff;
      milpFeasible = feasible;
    }
  }
  int numberPasses = 0;
  bool foundSolution = 0;
  while(isInteger && feasible ) {
    numberPasses++;

    //eventually give some information to user
    double time = CoinCpuTime();
    if(time - lastPeriodicLog > logFrequency_) {
      double lb = (model == NULL) ?si_->getObjValue():model->getBestPossibleObjValue();
      handler_->message(PERIODIC_MSG,messages_)
      <<time - timeBegin_<<cutoff
      <<lb
      <<CoinMessageEol;
      lastPeriodicLog = CoinCpuTime();
    }


    //setup the nlp
    int numberCutsBefore = cs.sizeRowCuts();

    bool fixed = true;

#ifdef OA_DEBUG
  std::cout<<"FX_P   : ";
#endif
  //Fix the variable which have to be fixed, after having saved the bounds
  for(int i = 0; i < numcols ; i++) {
      if(nlp_->isInteger(i)) {
        double value =  (model == NULL) ? si_->getColSolution()[i] : model->bestSolution()[i];
        if(saveColUb[i] < saveColLb[i] + cbcIntegerTolerance_)
          fixed = false;
      value = floor(value+0.5);
        value = max(saveColLb[i],value);
        value = min(value, saveColUb[i]);
        if(fabs(value) > 1e10) { std::cerr<<"FATAL ERROR: Variable taking a big value ("<<value
                                 <<") at optimium of LP relaxation, can not construct outer approximation. You should try running the problem with B-BB"<<std::endl;
         throw -1;
        }
#ifdef OA_DEBUG
        //         printf("xx %d at %g (bounds %g, %g)",i,value,nlp_->getColLower()[i],
        //                nlp_->getColUpper()[i]);
        std::cout<<(int)value;
#endif
        nlp_->setColLower(i,value);
        nlp_->setColUpper(i,value);
      }
    }
#ifdef OA_DEBUG
    std::cout<<std::endl;
#endif
    //Now solve the NLP get the cuts, and intall them in the local LP

    //  nlp_->turnOnIpoptOutput();
    nSolve_++;
    nlp_->resolve();
    if(nlp_->isProvenOptimal()) {
      nlp_->getOuterApproximation(cs);
      handler_->message(FEASIBLE_NLP, messages_)
      <<nlp_->getIterationCount()
      <<nlp_->getObjValue()<<CoinMessageEol;


#ifdef OA_DEBUG

      const double * colsol2 = nlp_->getColSolution();
      for(int i = 0 ; i < numcols; i++) {
        if(nlp_->isInteger(i)) {
          if(fabs(colsol2[i] - floor(colsol2[i] + 0.5) ) >
              cbcIntegerTolerance_) {
            std::cerr<<"Integer infeasible point (should not be), integer infeasibility for variable "<<i
            <<" is, "<<fabs(colsol2[i] - floor(colsol2[i] + 0.5))<<std::endl;
            //		    integerFeasible = 0;
          }
        }
      }
#endif

      if((nlp_->getObjValue() < cutoff) ) {
        handler_->message(UPDATE_UB, messages_)
        <<nlp_->getObjValue()
        <<CoinCpuTime()-timeBegin_
        <<CoinMessageEol;

        foundSolution = 1;
#if 1
        // Also store into solver
        if(babInfo) {
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
#endif
        // Update the cutoff
        cutoff = nlp_->getObjValue() *(1 - cbcCutoffIncrement_);
        // Update the lp solver cutoff
        si_->setDblParam(OsiDualObjectiveLimit, cutoff);
      }
    }
    else if(nlp_->isAbandoned() || nlp_->isIterationLimitReached()) {
      std::cerr<<"Unsolved NLP... exit"<<std::endl;
    }
    else {
      handler_->message(INFEASIBLE_NLP, messages_)
      <<nlp_->getIterationCount()
      <<CoinMessageEol;
      if(solveAuxiliaryProblem_)  //solve feasibility NLP and add corresponding outer approximation constraints
      {
        double *x = new double[numcols];
        int * ind = new int[numcols];
        int nInd = 0;
        for(int i = 0; i < numcols ; i++)
        {
          if(nlp_->isInteger(i)) {
            ind[nInd] = i;
            x[nInd++] = (model == NULL) ? si_->getColSolution()[i] : model->bestSolution()[i];
            //reset the bounds
            nlp_->setColBounds(i,saveColLb[i],saveColUb[i]);
          }
        }
        nlp_->getFeasibilityOuterApproximation(nInd, x, ind,cs);
        if(nlp_->isProvenOptimal())
        {
          assert((nlp_->getObjValue() > cbcIntegerTolerance_) );
          std::cout<<"Solved auxiliary infeasibility problem"<<std::endl;
        }
        else
        {
          std::cout<<"Failed to solve auxiliary infeasibility problem"<<std::endl;
          bool * s = const_cast<bool *> (&solveAuxiliaryProblem_);
          *s =0;
        }
        delete []x;
        delete[]ind;
      }
      else
        nlp_->getOuterApproximation(cs);
    }

    //install the cuts in si_ reoptimize and check for integer feasibility
    // Code taken and adapted from CbcModel::solveWithCuts
    int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
    if (numberCuts > 0) {
      CoinWarmStartBasis * basis
      = dynamic_cast<CoinWarmStartBasis*>(si_->getWarmStart()) ;
      assert(basis != NULL); // make sure not volume
      basis->resize(si_->getNumRows()+numberCuts,numcols + 1) ;
      for (int i = 0 ; i < numberCuts ; i++) {
        basis->setArtifStatus(si_->getNumRows()+i,
            CoinWarmStartBasis::basic) ;
      }

      const OsiRowCut ** addCuts = new const OsiRowCut * [numberCuts] ;
      for (int i = 0 ; i < numberCuts ; i++) {
        addCuts[i] = &cs.rowCut(i + numberCutsBefore) ;
      }
      si_->applyRowCuts(numberCuts,addCuts) ;

      delete [] addCuts ;
      if (si_->setWarmStart(basis) == false) {
        delete basis;
        throw CoinError("Fail setWarmStart() after cut installation.",
            "generateCuts","OACutGenerator2") ;
      }
      delete basis;
      si_->resolve();
      double objvalue = si_->getObjValue();
      //milpBound = max(milpBound, si_->getObjValue());
      feasible = (si_->isProvenOptimal() &&
          !si_->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
      //if value of integers are unchanged then we have to get out
      bool changed = //!nlp_->isProvenOptimal()//if nlp_ is infeasible solution has necessarily changed
          !feasible;//if lp is infeasible we don't have to check anything
      const double * nlpSol = nlp_->getColSolution();
      const double * lpSol = si_->getColSolution();
      for(int i = 0; i < numcols && !changed; i++) {
        if(nlp_->isInteger(i) && fabs(nlpSol[i] - lpSol[i])>0.001)
          changed = 1;
      }
      if(changed) {
        isInteger = 1;
        for(int i = 0 ; i < numcols ; i++) {
          if(nlp_->isInteger(i)) {
            if(fabs(lpSol[i] - floor(lpSol[i] + 0.5) ) > cbcIntegerTolerance_ ) {
              isInteger = 0;
              break;
            }
          }
        }
      }
      else {
        isInteger = 0;
        //	  if(!fixed)//fathom on bounds
        milpBound = 1e200;
      }
#ifdef OA_DEBUG

      printf("Obj value after cuts %g %d rows\n",si_->getObjValue(),
          numberCuts) ;
#endif
      //do we perform a new local search ?
      if(milpOptimal && feasible && !isInteger &&
          nLocalSearch_ < maxLocalSearch_ &&
          numberPasses < maxLocalSearchPerNode_ &&
          localSearchNodeLimit_ > 0 &&
          CoinCpuTime() - timeBegin_ < maxLocalSearchTime_) {
         std::cout<<"Perform new local search"<<std::endl;
        if(clp==NULL) {
          nLocalSearch_++;

          setTimeLimit(maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());
          setCutoff(cutoff);
          si_->branchAndBound();
          //Did we find any integer feasible solution
          milpBound = siBestObj();
          //If integer solution is the same as nlp solution problem is solved
          changed = !nlp_->isProvenOptimal();
          for(int i = 0; i < numcols && !changed; i++) {
            if(nlp_->isInteger(i) && fabs(nlp_->getColSolution()[i] - si_->getColSolution()[i])>cbcIntegerTolerance_)
              changed = 1;
          }
          if(changed) {
            feasible = (milpBound < cutoff);
            std::cout<<"milp bound "<<milpBound<<" cutoff "<<cutoff<<std::endl;
          }
          else
           {
            feasible = 0;
            milpBound = 1e50;
           }
          isInteger = si_->getFractionalIndices().size() == 0;
        }/** endif solved by branchAndBound capable OsiInterface*/
        else {
          if(model)
	    {
	      OAModel = NULL;	    
	      delete model;
	    }
          model = new CbcModel(*clp); // which clones
	  OAModel = model;
          OsiBabSolver empty;
          model->solver()->setAuxiliaryInfo(&empty);
          //Change Cbc messages prefixes
          strcpy(model->messagesPointer()->source_,"OaCbc");
          model->solver()->messageHandler()->setLogLevel(0);

          if(!strategy)
            strategy = new CbcStrategyDefault(1,0,0, subMilpLogLevel_);
          model->setLogLevel(subMilpLogLevel_);
          model->setStrategy(*strategy);
          model->setMaximumNodes(localSearchNodeLimit_);
          model->setMaximumSeconds(maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());
          model->setCutoff(cutoff);
          model->branchAndBound();
          nLocalSearch_++;
          milpBound = siBestObj(model);
          handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<model->getNodeCount()<<model->getIterationCount()<<CoinMessageEol;
          feasible =  (milpBound < cutoff);
          isInteger = model->getSolutionCount();
          if(model->isProvenOptimal() || model->isProvenInfeasible()) {
            bool changed = 0;		  //If integer solution is the same as nlp solution problem is solved
            for(int i = 0; feasible && isInteger && i < numcols && !changed; i++) {
              if(nlp_->isInteger(i) && fabs(nlp_->getColSolution()[i] - model->bestSolution()[i])>cbcIntegerTolerance_)
                changed = 1;
            }
          if(!changed) {
            feasible = 0;
            milpBound = 1e50;
           }
            milpFeasible = feasible;
          }
          else {
            std::cout<<"Exiting on time limit"<<std::endl;
            milpOptimal = 0;
            //feasible = 1;
            break;
          }

        }/** endif solved by clp/cbc*/

        if(milpBound < cutoff)
          handler_->message(UPDATE_LB, messages_)
          <<milpBound<<CoinCpuTime() - timeBegin_
          <<CoinMessageEol;
        else
        {
          milpBound = 1e50;
          feasible = 0;
          handler_->message(OASUCCESS, messages_)
          <<CoinCpuTime() - timeBegin_ <<CoinMessageEol;
        }
      }/** endif localSearch*/
      else if(model!=NULL)/** have to delete model for next iteration*/
      { 
	delete model;
	OAModel = NULL;
        model=NULL;
      }

    }



  }
#ifdef OA_DEBUG
    std::cout<<"-------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"OA cut generation finished"<<std::endl;
    std::cout<<"Generated "<<cs.sizeRowCuts()<<std::endl;
    if(foundSolution)
      std::cout <<"Found NLP-integer feasible solution of  value : "<<cutoff<<std::endl;
    std::cout<<"Current MILP lower bound is : "<<milpBound<<std::endl;
    std::cout<<"-------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"Stopped because : isInteger "<<isInteger<<", feasible "<<feasible<<std::endl<<std::endl;
#endif
    //Transmit the bound found by the milp
    {
      if(milpBound>-1e100)
      {
#if 1
        // Also store into solver
        if (babInfo)
          babInfo->setMipBound(milpBound);
#endif

      }
    }  //Clean everything :

  //free model
  if(model!= NULL) {
    delete model;
  }



  //  Reset bounds in the NLP

  for(int i = 0; i < numcols ; i++) {
    if(nlp_->isInteger(i)) {
      nlp_->setColBounds(i,saveColLb[i],saveColUb[i]);
    }
  }
#ifndef NO_NULL_SI
  if(deleteSi) {
    delete si_;
    si_ = NULL;
  }
  else
#endif 
   {
    if(leaveSiUnchanged_) {
      int nRowsToDelete = si_->getNumRows() - originalRowNumber;
      int * rowsToDelete = new int[nRowsToDelete];
      for(int i = 0 ; i < nRowsToDelete ; i++) {
        rowsToDelete[i] = i + originalRowNumber;
      }
      si_->deleteRows(nRowsToDelete, rowsToDelete);
      delete [] rowsToDelete;
      if(si_==&si) {
        si_->setDblParam(OsiDualObjectiveLimit, saveCutoff);
        if(si_->setWarmStart(saveWarmStart)==false) {
          throw CoinError("Fail restoring the warm start at the end of procedure",
              "generateCuts","OACutGenerator2") ;
        }
        delete saveWarmStart;
      }
    }
}

delete [] saveColLb;
delete [] saveColUb;
if(strategy && ! strategy_)
  delete strategy;
}
}
