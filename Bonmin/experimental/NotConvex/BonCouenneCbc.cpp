// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006


//Bonmin header files
#include "BonminConfig.h"
#include "BonCbc.hpp"
#include "BonCbcLpStrategy.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonOsiTMINLPInterface.hpp"

#include "BonCbcParam.hpp"

// AW
#include "BonCurvBranching.hpp"
#include "BonQPStrongBranching.hpp"
#include "BonLpStrongBranching.hpp"

//OA machinery
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaFeasChecker.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonEcpCuts.hpp"

// Cbc Header file
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"

//Osi Header files
#include "OsiClpSolverInterface.hpp"
#include "OsiAuxInfo.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif



// Cut generators
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"

// Node selection
#include "CbcCompareUser.hpp"
#include "CbcCompareActual.hpp"

#include "CbcBranchUser.hpp"


// Code to enable user interuption
static CbcModel * currentBranchModel = NULL; //pointer to the main b&b
static Bonmin::OACutGenerator2 * currentOA = NULL; //pointer to the OA generator
CbcModel * OAModel = NULL; // pointer to the submip if using Cbc
bool BonminAbortAll = false;
//#ifdef COIN_HAS_CPX
//OsiCpxSolverInterface * CpxModel = NULL;//pointer to the submip if using cplex
//#endif
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
  /** Constructor.*/
  CouenneBab::CouenneBab():
    Bab()
  {}

  /** Destructor.*/
  CouenneBab::~CouenneBab()
  {
  }
  /** Perform a branch-and-bound on given IpoptInterface using passed parameters.*/
  void
  CouenneBab::branchAndBound(OsiTMINLPInterface *nlpSolver,
      const BonminCbcParam &par)
  {

    //Now set-up b&b
    OsiSolverInterface * si;


    nlpSolver->messageHandler()->setLogLevel(par.nlpLogLevel);

    if (par.algo >= 0 && par.algo <= 3) //Do something classic
    {
      Bab::branchAndBound(nlpSolver, par);
      return;
    }

    if(!dynamic_cast<BonCouenneInterface *>(nlpSolver))
      {
	std::cerr<<"Can not do Couenne algorithm without a Couenne interface."
		 <<std::endl;
      }

    si = new OsiClpSolverInterface;
    nlpSolver->extractLinearRelaxation(*si);
    // say bound dubious, does cuts at solution
    OsiBabSolver * extraStuff = new OsiBabSolver(0);
    si->setAuxiliaryInfo(extraStuff);
    delete extraStuff;


    CbcModel model(*si);


    //Setup likely milp cut generator
    CglGomory miGGen;

    CglProbing probGen;
    probGen.setUsingObjective(true);
    probGen.setMaxPass(3);
    probGen.setMaxProbe(100);
    probGen.setMaxLook(50);

    CglKnapsackCover knapsackGen;
    CglMixedIntegerRounding mixedGen;

    //Setup OA generators

    EcpCuts ecpGen(nlpSolver);
    ecpGen.parameter().global_ = par.oaCutsGlobal;
    ecpGen.parameter().addOnlyViolated_ = par.addOnlyViolatedOa;


    int numGen = 0;

    if(par.filmintCutsFrequency != 0){
      model.addCutGenerator(&ecpGen,par.filmintCutsFrequency,"Filmint cutting planes");
      numGen++;
    }

    if (par.migFreq != 0) {
      model.addCutGenerator(&miGGen,par.migFreq,"GMI");
      numGen++;
    }
      if (par.probFreq != 0) {
        model.addCutGenerator(&probGen,par.probFreq,"Probing");
        numGen++;
      }
      if (par.coverFreq != 0) {
        model.addCutGenerator(&knapsackGen,par.coverFreq,"covers");
        numGen++;
      }
      if (par.mirFreq != 0) {
        model.addCutGenerator(&mixedGen,par.mirFreq,"MIR");
        numGen++;
      }

    //Set true branch-and-bound parameters
    model.messageHandler()->setLogLevel(par.bbLogLevel);
    if (par.algo > 0)
      model.solver()->messageHandler()->setLogLevel(par.lpLogLevel);


    //   model.setMaxFailure(par.maxFailures);
    //   model.setMaxInfeasible(par.maxInfeasible);

    //TODO add the couenne objects to model


    //Pass over user set branching priorities to Cbc
    {
      //set priorities, prefered directions...
      const int * priorities = nlpSolver->getPriorities();
      const double * upPsCosts = nlpSolver->getUpPsCosts();
      const double * downPsCosts = nlpSolver->getDownPsCosts();
      const int * directions = nlpSolver->getBranchingDirections();
      bool hasPseudo = (upPsCosts!=NULL);
      model.findIntegers(true,hasPseudo);
      OsiObject ** simpleIntegerObjects = model.objects();
      int numberObjects = model.numberObjects();
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

    // Now pass user set Sos constraints (code inspired from CoinSolve.cpp)
    const TMINLP::SosInfo * sos = nlpSolver->model()->sosConstraints();
    if (!par.disableSos && sos && sos->num > 0) //we have some sos constraints
    {
      const int & numSos = sos->num;
      CbcObject ** objects = new CbcObject*[numSos];
      const int * starts = sos->starts;
      const int * indices = sos->indices;
      const char * types = sos->types;
      const double * weights = sos->weights;
      //verify if model has user set priorities
      bool hasPriorities = false;
      const int * varPriorities = nlpSolver->getPriorities();
      int numberObjects = model.numberObjects();
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
        objects[i] = new CbcSOS(&model, length, &indices[start],
            &weights[start], i, types[i]);

        objects[i]->setPriority(10);
        if (hasPriorities && sosPriorities && sosPriorities[i]) {
          objects[i]->setPriority(sosPriorities[i]);
        }
      }
      model.addObjects(numSos, objects);
      for (int i = 0 ; i < numSos ; i++)
        delete objects[i];
      delete [] objects;
    }


    replaceIntegers(model.objects(), model.numberObjects());

    model.setPrintFrequency(par.logInterval);

    model.setDblParam(CbcModel::CbcCutoffIncrement, par.cutoffDecr);

    model.setCutoff(par.cutoff);
    //  model.setBestObjectiveValue(par.cutoff);

    model.setDblParam(CbcModel::CbcAllowableGap, par.allowableGap);
    model.setDblParam(CbcModel::CbcAllowableFractionGap, par.allowableFractionGap);

    // Definition of node selection strategy
    CbcCompareObjective compare0;
    CbcCompareDepth compare1;
    CbcCompareUser compare2;
    if (par.nodeSelection==0) {
      model.setNodeComparison(compare0);
    }
    else if (par.nodeSelection==1) {
      model.setNodeComparison(compare1);
    }
    else if (par.nodeSelection==2) {
      compare2.setWeight(0.0);
      model.setNodeComparison(compare2);
    }
    else if (par.nodeSelection==3) {
      model.setNodeComparison(compare2);
    }

    model.setNumberStrong(par.numberStrong);

    model.setNumberBeforeTrust(par.minReliability);

    model.setNumberPenalties(8);

    model.setDblParam(CbcModel::CbcMaximumSeconds, par.maxTime);

    model.setMaximumNodes(par.maxNodes);

    model.setIntegerTolerance(par.intTol);


    // Redundant definition of default branching (as Default == User)
    CbcBranchUserDecision branch;


    if(par.varSelection == OsiTMINLPInterface::CURVATURE_ESTIMATOR){
    // AW: Try to set new chooseVariable object
      BonCurvBranching chooseVariable(nlpSolver);
      chooseVariable.setNumberStrong(model.numberStrong());
      branch.setChooseMethod(chooseVariable);
    }
    else if(par.varSelection == OsiTMINLPInterface::QP_STRONG_BRANCHING){
      model.solver()->findIntegersAndSOS(false);
      BonQPStrongBranching chooseVariable(nlpSolver);
      chooseVariable.setNumberStrong(model.numberStrong());
      branch.setChooseMethod(chooseVariable);
    }
    else if(par.varSelection == OsiTMINLPInterface::LP_STRONG_BRANCHING){
      model.solver()->findIntegersAndSOS(false);
      LpStrongBranching chooseVariable(nlpSolver);
      chooseVariable.setMaxCuttingPlaneIter(par.numEcpRounds);
      chooseVariable.setNumberStrong(model.numberStrong());
      branch.setChooseMethod(chooseVariable);
    }
    else if(par.varSelection == OsiTMINLPInterface::NLP_STRONG_BRANCHING){
      const bool solve_nlp = true;
      model.solver()->findIntegersAndSOS(false);
      BonQPStrongBranching chooseVariable(nlpSolver, solve_nlp);
      chooseVariable.setNumberStrong(model.numberStrong());
      branch.setChooseMethod(chooseVariable);
    }
    else if(par.varSelection == OsiTMINLPInterface::OSI_SIMPLE){
      model.solver()->findIntegersAndSOS(false);
      OsiChooseVariable choose(model.solver());
      branch.setChooseMethod(choose);
    }
    else if(par.varSelection == OsiTMINLPInterface::OSI_STRONG){
      model.solver()->findIntegersAndSOS(false);
      OsiChooseStrong choose(model.solver());
      choose.setNumberBeforeTrusted(par.minReliability);
      choose.setNumberStrong(par.numberStrong);
      branch.setChooseMethod(choose);
    }
    model.setBranchingMethod(&branch);

    //Get the time and start.
    model.initialSolve();

    continuousRelaxation_ =model.solver()->getObjValue();
    if (par.algo == 0)//Set warm start point for Ipopt
    {
      const double * colsol = model.solver()->getColSolution();
      const double * duals = model.solver()->getRowPrice();
      model.solver()->setColSolution(colsol);
      model.solver()->setRowPrice(duals);
    }

#ifdef SIGNAL
    CoinSighandler_t saveSignal=SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
#endif

    currentBranchModel = &model;

    model.branchAndBound();
    numNodes_ = model.getNodeCount();
    bestObj_ = model.getObjValue();
    bestBound_ = model.getBestPossibleObjValue();
    mipIterationCount_ = model.getIterationCount();

    bool hasFailed = false;
    if (par.algo==0)//Did we continue branching on a failure
    {
      CbcNlpStrategy * nlpStrategy = dynamic_cast<CbcNlpStrategy *>(model.strategy());
      if (nlpStrategy)
        hasFailed = nlpStrategy->hasFailed();
      else
        throw -1;
    }
    else
      hasFailed = nlpSolver->hasContinuedOnAFailure();


    if (hasFailed) {
      std::cout<<"************************************************************"<<std::endl
      <<"WARNING : Optimization failed on an NLP during optimization"
      <<"\n (no optimal value found within tolerances)."<<std::endl
      <<"Optimization was not stopped because option \n"
      <<"\"nlp_failure_behavior\" has been set to fathom but"
      <<" beware that reported solution may not be optimal"<<std::endl
      <<"************************************************************"<<std::endl;
    }

    if (model.bestSolution()) {
      if (bestSolution_)
        delete [] bestSolution_;
      bestSolution_ = new double[nlpSolver->getNumCols()];
      CoinCopyN(model.bestSolution(), nlpSolver->getNumCols(), bestSolution_);
    }
    if (!model.status()) {
      if (bestSolution_)
        mipStatus_ = FeasibleOptimal;
      else
        mipStatus_ = ProvenInfeasible;
    }
    else {
      if (bestSolution_)
        mipStatus_ = Feasible;
      else
        mipStatus_ = NoSolutionKnown;
    }
    delete si;
    std::cout<<"Finished"<<std::endl;
    if (strategy)
      delete strategy;

  }

}
