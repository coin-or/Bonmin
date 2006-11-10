// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006


//Bonmin heaader files
#include "BonminConfig.h"
#include "BonCbc.hpp"
#include "BonCbcLpStrategy.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonOsiTMINLPInterface.hpp"

#include "BonCbcParam.hpp"

// AW
#include "BonChooseVariable.hpp"

//OA machinery
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOACutGenerator.hpp"


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
//#ifdef COIN_HAS_CPX
//OsiCpxSolverInterface * CpxModel = NULL;//pointer to the submip if using cplex
//#endif

#include "CoinSignal.hpp"

extern "C"
{
  static void signal_handler(int whichSignal) {
    if (currentBranchModel!=NULL)
      currentBranchModel->setMaximumNodes(0); // stop at next node
    if (OAModel!=NULL)
      OAModel->setMaximumNodes(0); // stop at next node
    if (currentOA!=NULL)
      currentOA->setMaxLocalSearchTime(0.); // stop OA
    return;
  }
}

namespace Bonmin
{
  /** Constructor.*/
  Bab::Bab():
      bestSolution_(NULL),
      mipStatus_(),
      bestObj_(1e200),
      bestBound_(-1e200),
      continuousRelaxation_(-DBL_MAX),
      numNodes_(0),
      mipIterationCount_(0)
  {}

  /** Destructor.*/
  Bab::~Bab()
  {
    if (bestSolution_) delete [] bestSolution_;
    bestSolution_ = NULL;
  }
  /** Perform a branch-and-bound on given IpoptInterface using passed parameters.*/
  void
  Bab::branchAndBound(OsiTMINLPInterface *nlpSolver,
      const BonminCbcParam &par)
  {

    //Now set-up b&b
    OsiSolverInterface * si;


    nlpSolver->messageHandler()->setLogLevel(par.nlpLogLevel);

    if (par.algo > 0) //OA based
    {
      si = new OsiClpSolverInterface;
      nlpSolver->extractLinearRelaxation(*si);
      // say bound dubious, does cuts at solution
      OsiBabSolver * extraStuff = new OsiBabSolver(3);
      si->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
    else {
      si = nlpSolver;
      nlpSolver->ignoreFailures();
      OsiBabSolver * extraStuff = new OsiBabSolver(2);
      si->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
    CbcModel model(*si);



    if (par.algo==0)//Switch off some feasibility checks and switch on specific nodes info
    {
      int specOpt = model.specialOptions();
      specOpt = 16;
      model.setSpecialOptions(specOpt);
      CbcNlpStrategy strat(par.maxFailures, par.maxInfeasible, par.failureBehavior);
      model.setStrategy(strat);
    }

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

    //Resolution of nlp relaxations
    OACutGenerator oaGen(nlpSolver);
    oaGen.setMaxDepth(100000);
    oaGen.setLogLevel(par.oaLogLevel);


    //Outer approximation iterations
    OsiSolverInterface * localSearchSolver=NULL;
    if (par.milpSubSolver <= 1)/* use cbc */
    {
      localSearchSolver = model.solver();
      //localSearchSolver->messageHandler()->setLogLevel(0);
    }
    else if (par.milpSubSolver ==2) /* try to use cplex */
    {
#ifdef COIN_HAS_CPX
      OsiCpxSolverInterface * cpxSolver = new OsiCpxSolverInterface;
      localSearchSolver = cpxSolver;
      nlpSolver->extractLinearRelaxation(*localSearchSolver);
#else

      std::cerr<<"You have set an option to use CPLEX as the milp subsolver in oa decomposition."<<std::endl
      <<"but apparently CPLEX is not configured to be used in bonmin, see the manual for configuring CPLEX"<<std::endl;
      throw -1;
#endif

    }
    CbcStrategy * strategy = NULL;
    if (par.milpSubSolver == 1) {
      strategy = new CbcOaStrategy(par.milpSubSolver_migFreq,
          par.milpSubSolver_probFreq,
          par.milpSubSolver_mirFreq,
          par.milpSubSolver_coverFreq,
          par.milpSubSolver_minReliability,
          par.milpSubSolver_numberStrong,
          par.milpSubSolver_nodeSelection,
          par.intTol,
          par.milpLogLevel
                                  );
    }
    OACutGenerator2 oaDec(nlpSolver, localSearchSolver, strategy, par.cutoffDecr, par.intTol, 0,1);
    if (par.algo>0) {
      oaDec.setLocalSearchNodeLimit(1000000);
      oaDec.setMaxLocalSearch(100000);
      oaDec.setMaxLocalSearchPerNode(10000);
      oaDec. setMaxLocalSearchTime(min(par.maxTime,par.oaDecMaxTime));
      oaDec.setLogLevel(par.oaLogLevel);
      oaDec.setLogFrequency(par.oaLogFrequency);
      oaDec.setSubMilpLogLevel(par.milpLogLevel);
    }
    //Setup solver for checking validity of integral solutions
    OACutGenerator2 feasCheck(nlpSolver, model.solver(),
        NULL,
        par.cutoffDecr, par.intTol,
        0, 0);
    if (par.algo>0) {
      feasCheck.setLocalSearchNodeLimit(0);
      feasCheck.setMaxLocalSearch(0);
      feasCheck.setMaxLocalSearchPerNode(100000);
    }
    DummyHeuristic oaHeu(model, nlpSolver);

    if (par.algo>0) {
      int numGen = 0;
      if (par.nlpSolveFrequency != 0) {
        model.addCutGenerator(&oaGen,par.nlpSolveFrequency,"Outer Approximation Supporting Hyperplanes for NLP optimum");
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
      if (par.oaDecMaxTime>0.) {
        model.addCutGenerator(&oaDec,1,"Outer Approximation local enumerator");
        OACutGenerator2 * oaDecCopy = dynamic_cast<OACutGenerator2 *>
            (model.cutGenerators()[numGen]->generator());
        assert(oaDecCopy);
        currentOA = oaDecCopy;
        numGen++;
      }
      model.addCutGenerator(&feasCheck,1,"Outer Approximation feasibility checker",false,true);
      numGen++;

      model.addHeuristic(&oaHeu);
    }

    //Set true branch-and-bound parameters
    model.messageHandler()->setLogLevel(par.bbLogLevel);
    if (par.algo > 0)
      model.solver()->messageHandler()->setLogLevel(par.lpLogLevel);


    //   model.setMaxFailure(par.maxFailures);
    //   model.setMaxInfeasible(par.maxInfeasible);

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


    if(par.varSelection == 3){
    // AW: Try to set new chooseVariable object
      BonChooseVariable chooseVariable(nlpSolver);
      chooseVariable.setNumberStrong(model.numberStrong());
      branch.setChooseMethod(chooseVariable);}
     else if(par.varSelection == 4){
      OsiChooseVariable choose(model.solver());
      branch.setChooseMethod(choose);
    }
    else if(par.varSelection == 5){
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

    CoinSighandler_t saveSignal=SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);


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


    if (par.algo > 0)
      delete si;
#ifdef COIN_HAS_CPX

    if (par.milpSubSolver ==1)
      delete localSearchSolver;
#endif
    std::cout<<"Finished"<<std::endl;
    if (strategy)
      delete strategy;

  }


  /** return the best known lower bound on the objective value*/
  double
  Bab::bestBound()
  {
    if (mipStatus_ == FeasibleOptimal) return bestObj_;
    else if (mipStatus_ == ProvenInfeasible) return 1e200;
    else return bestBound_;
  }
}
