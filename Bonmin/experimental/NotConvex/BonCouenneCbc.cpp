// (C) Copyright International Business Machines Corporation (IBM) and Carnegie Mellon University 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belloti, Carnegie Mellon University
// Pierre Bonami, International Business Machines Corporation

//
// Date : 02/14/2007


//Couenne Bonmin interface
#include "BonCouenneCbc.hpp"
#include "BonCouenneInterface.hpp"
#include "BonNlpHeuristic.hpp"
//#include "BonCouenneConvexCuts.hpp"

//Couenne
#include "CouenneObject.hpp"
#include "CouenneChooseVariable.hpp"

//Bonmin header files
#include "BonminConfig.h"
#include "BonCbcLpStrategy.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonOsiTMINLPInterface.hpp"

#include "BonCbcParam.hpp"

// AW
#include "BonCurvBranching.hpp"
#include "BonQPStrongBranching.hpp"
#include "BonLpStrongBranching.hpp"



// Cbc Header file
#include "CbcModel.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"

//Osi Header files
#include "OsiClpSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiBranchingObject.hpp"
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

bool checkNLP (CglCutGenerator *, const double *, double &);


// Code to enable user interuption
static CbcModel * currentBranchModel = NULL; //pointer to the main b&b
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
      std::cerr<<"User forced interruption"<<std::endl;
      exit(0);
    }
    if (currentBranchModel!=NULL)
      currentBranchModel->setMaximumNodes(0); // stop at next node
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
    CouenneInterface* ci = dynamic_cast<CouenneInterface *>(nlpSolver);

    nlpSolver->messageHandler()->setLogLevel(par.nlpLogLevel);

    if (par.algo >= OsiTMINLPInterface::B_BB && 
	par.algo <= OsiTMINLPInterface::B_Hyb) //Do something classic
    {
      Bab::branchAndBound (nlpSolver, par);
      return;
    }

    if (ci == NULL)
	std::cerr << "Can not run Couenne without a Couenne interface."
		  << std::endl;

    si = new OsiClpSolverInterface;

    nlpSolver -> extractLinearRelaxation (*si);

    // [PBelotti] Change solver type to force NLP-feasibility check
    // rather than simple integrality
    //    OsiBabSolver * extraStuff = new OsiBabSolver(3);

    // say bound dubious, does cuts at solution
    OsiBabSolver * extraStuff = new OsiBabSolver(0);

    // as per instructions by John Forrest, to get changed bounds
    extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);

    si -> setAuxiliaryInfo (extraStuff);
    delete extraStuff;

    //TODO : Switch to OsiObejcts
    si->findIntegersAndSOS(false);
    int numberIntegerObjects = si->numberObjects() > 0;
    {
      const CouenneProblem * couenneProb = ci->couenneProb();
      int numAuxs = couenneProb->nAuxs();
      OsiObject ** objects = new OsiObject*[numAuxs];
      int nobj = 0;

      for (int i = 0 ; i < numAuxs; i++) // for each aux variable

	// if this variable is associated with a nonlinear function
	//	if (couenneProb -> Aux (i) -> Image () -> Linearity () > LINEAR) 
	  {
	    /*
	    printf ("creating CouenneObject for ");

	    couenneProb -> Aux (i) ->             print (std::cout); printf (" := ");
	    couenneProb -> Aux (i) -> Image () -> print (std::cout); printf ("\n");
	    */
	    // then we may have to branch on it
	    objects [nobj] = new CouenneObject (couenneProb -> Aux (i));
	    objects [nobj++] -> setPriority (1);
	  }

      si -> addObjects (nobj, objects);
//       for(int i = 0 ; i < nobj ; i++){
// 	delete objects[i];
//       }
//       delete objects;
    }
    

    // create linear model based on the initial linear convexification

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

    //Setup Convexifier generators

    CouenneCutGenerator *ecpGen = ci -> couenneCg ();

    ecpGen -> setBabPtr (this);


    int numGen = 0;

    if (par.couenneCutsFrequency != 0) {
      model.addCutGenerator (ecpGen, par.couenneCutsFrequency, "Couenne cutting planes");
      numGen++;
    }
    /*
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
    */

    /*Setup heuristic to solve nlp problems.*/
    Ipopt::SmartPtr<Ipopt::OptionsList> Options = nlpSolver->retrieve_options();
    int doNlpHeurisitic = 0;
    Options->GetEnumValue("nlp_heuristic_solves", doNlpHeurisitic, "couenne.");
    if(doNlpHeurisitic)
    {
      NlpSolveHeuristic nlpHeuristic(model, *nlpSolver, false);
      nlpHeuristic.setMaxNlpInf(1e10);
      model.addHeuristic(&nlpHeuristic);
      //Set true branch-and-bound parameters
      model.messageHandler()->setLogLevel(par.bbLogLevel);
    }
    if (par.algo > 0)
      model.solver()->messageHandler()->setLogLevel(par.lpLogLevel);


    //   model.setMaxFailure(par.maxFailures);
    //   model.setMaxInfeasible(par.maxInfeasible);


    //Pass over user set branching priorities to Cbc
    if(0)
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

    CouenneChooseVariable choose (model.solver(), 
				  const_cast<CouenneProblem *> (ci -> couenneProb ()));
    branch.setChooseMethod(choose);

    model.setBranchingMethod(&branch);

    //Get the time and start.
    model.initialSolve();

    continuousRelaxation_ = model.solver() -> getObjValue ();
    double * nlpSol = NULL;
    if(nlpSolver->isProvenOptimal()){
      //Look at the solution of nlp if it satisfies all integer branching object (if any) set it as best solution.
      
      nlpSol = new double[model.solver()->getNumCols()];
      CoinCopyN(nlpSolver->getColSolution(), nlpSolver->getNumCols(), nlpSol);
      CouenneInterface * couenne = dynamic_cast<CouenneInterface *>
        (nlpSolver);
      if(couenne){
        couenne->couenneProb()->getAuxs(nlpSol);
      }
      
      OsiObject ** objects = si->objects();
      OsiBranchingInformation info(model.solver(), true, false);
      info.solution_ = nlpSol;
      bool feasible = true;
      for(int i = 0 ; i < numberIntegerObjects ; i++){
        int dummy;
        if(objects[i]->infeasibility(&info,dummy) > par.intTol)
        {
          feasible = false;
          break;
        }
      }
      if(feasible){
        double nlpObj = nlpSolver->getObjValue();
        model.setSolutionCount(1);
        model.setObjValue(nlpObj);
        model.setCutoff(nlpObj);
      }
      else
        delete [] nlpSol;
    }
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

    //////////////////////////////////////////////////////////

    model.branchAndBound ();

    //////////////////////////////////////////////////////////

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

      double obj = 0;

      if (bestSolution_)
	if (!checkNLP (model.cutGenerator (0) -> generator (), bestSolution_, obj))
	  printf ("### checkNLP: solution is wrong\n");

      if (bestSolution_)
        delete [] bestSolution_;

      bestSolution_ = new double[nlpSolver->getNumCols()];
      CoinCopyN(model.bestSolution(), nlpSolver->getNumCols(), bestSolution_);

      //Check solution validity
      double violation = nlpSolver->getConstraintsViolation(bestSolution_, obj);

      if (fabs (violation) > 1e-5)
	std::cout << "Solution found violates non linear constraints by: "
		  << violation << std::endl << "Objective: " << obj
		  << " (was "<< model.getObjValue () << " returned by Cbc)." << std::endl;
      /*      
      printf ("lp solution: {");
      for (int i=0; i < ecpGen -> getnvars (); i++)
	printf ("%.3f ", bestSolution_ [i]);
      printf ("}\n");
      */
    }
    else if(nlpSol)
    {
      if (bestSolution_)
        delete [] bestSolution_;
      bestSolution_ = nlpSol;
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


    if (0) { // print some statistics in LaTeX format

      int nr, nt;
      double st;

      CouenneCutGenerator *cg = dynamic_cast <CouenneCutGenerator *> 
	(model.cutGenerator (0) -> generator ());

      cg -> getStats (nr, nt, st);

      printf ("::: %6d & %6d & %6d & %6d & %6d & %6d & %8.3f & ", 
	      cg -> Problem () -> nVars (), 
	      cg -> Problem () -> nIntVars(), 
	      cg -> Problem () -> nNLCons (),
	      cg -> Problem () -> nAuxs (),
	      nr, nt, st);
    }
  }
}
