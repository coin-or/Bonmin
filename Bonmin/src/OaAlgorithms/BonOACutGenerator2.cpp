// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 05/26/2005
//#define OA_DEBUG

#include "BonOACutGenerator2.hpp"
#include "BonminConfig.h"

#include "OsiClpSolverInterface.hpp"

#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#include "OsiAuxInfo.hpp"



namespace Bonmin
{


/// Default constructor
  OACutGenerator2::OACutGenerator2():
      CglCutGenerator(),
      OaDecompositionHelper()
   {
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
      OaDecompositionHelper(nlp,si,
                          strategy, cbcCutoffIncrement,
                          cbcIntegerTolerance, leaveSiUnchanged)
  {
  }

  OACutGenerator2::~OACutGenerator2()
  {
  }

  /// cut generation method
  void
  OACutGenerator2::generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
      const CglTreeInfo info) const
  {
    double lastPeriodicLog= CoinCpuTime();

    if (nlp_ == NULL) {
      std::cerr<<"Error in cut generator for outer approximation no NLP ipopt assigned"<<std::endl;
      throw -1;
    }

    // babInfo is used to communicate with the b-and-b solver (Cbc or Bcp).
    OsiBabSolver * babInfo = dynamic_cast<OsiBabSolver *> (si.getAuxiliaryInfo());

    const int numcols = nlp_->getNumCols();

    //Get the continuous solution
    const double *colsol = si.getColSolution();


    //Check integer infeasibility
    bool isInteger = integerFeasible(colsol, numcols);

    SubMipSolver * subMip = NULL;

    if (!isInteger) {
      if (nLocalSearch_<parameters_.maxLocalSearch_ &&
          parameters_.localSearchNodeLimit_ > 0 &&
          CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_)//do a local search
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
    solverManip nlpManip(nlp_, false, false, true, false);

    solverManip * lpManip = NULL; 
    OsiSolverInterface *lp;
    if(lp_ != NULL){
      if(lp_!=&si){
        lpManip = new solverManip(lp_, true, false, false, true);
        lpManip->cloneOther(si);
      }
      else{
#if 0
        throw CoinError("Not allowed to modify si in a cutGenerator",
          "OACutGenerator2","generateCuts");
#else
         lpManip = new solverManip(lp_, true, leaveSiUnchanged_, true, true);
#endif
      }
      lp = lp_;
    }
    else{
      lpManip = new solverManip(si);
      lp = lpManip->si();
    }

    bool milpOptimal = 1;


    double milpBound = -DBL_MAX;
    bool milpFeasible = 1;
    bool feasible = 1;

    if (subMip)//Perform a local search
    {
        subMip->performLocalSearch(cutoff, parameters_.subMilpLogLevel_, 
                                   (parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime()) /* time limit */,
                                   parameters_.localSearchNodeLimit_);   
        milpBound = subMip->lowBound();
        milpOptimal = subMip->optimal();

        feasible = milpBound < cutoff;
        milpFeasible = feasible;
        isInteger = subMip->getLastSolution() != NULL;
        nLocalSearch_++;

        if (milpOptimal)
          handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
        else
        {
          handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
        }
    }
    int numberPasses = 0;
    bool foundSolution = 0;
    while (isInteger && feasible ) {
      numberPasses++;

      //after a prescribed elapsed time give some information to user
      double time = CoinCpuTime();
      if (time - lastPeriodicLog > parameters_.logFrequency_) {
        double lb = (subMip == NULL) ?lp->getObjValue() : subMip->getLowerBound();
        handler_->message(PERIODIC_MSG,messages_)
        <<time - timeBegin_<<cutoff
        <<lb
        <<CoinMessageEol;
        lastPeriodicLog = CoinCpuTime();
      }


      //setup the nlp
      int numberCutsBefore = cs.sizeRowCuts();

    //Fix the variable which have to be fixed, after having saved the bounds
    colsol = (subMip == NULL) ? lp->getColSolution() : subMip->getLastSolution();
    nlpManip.fixIntegers(colsol);


    if(solveNlp(babInfo, cutoff)){
      //nlp solved and feasible
      // Update the cutoff
      cutoff = nlp_->getObjValue() *(1 - parameters_.cbcCutoffIncrement_);
      // Update the lp solver cutoff
      lp->setDblParam(OsiDualObjectiveLimit, cutoff);
    }

      
    // Get the cuts outer approximation at the current point
    nlp_->getOuterApproximation(cs);


    int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
    if (numberCuts > 0)
      lpManip->installCuts(cs, numberCuts);

    lp->resolve();
    double objvalue = lp->getObjValue();
    //milpBound = max(milpBound, lp->getObjValue());
    feasible = (lp->isProvenOptimal() &&
		!lp->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
    //if value of integers are unchanged then we have to get out
    bool changed = !feasible;//if lp is infeasible we don't have to check anything
    if(!changed)
      changed = nlpManip.isDifferentOnIntegers(lp->getColSolution());
    if (changed) {
      isInteger = integerFeasible(lp->getColSolution(), numcols);
    }
    else {
      isInteger = 0;
      //	  if(!fixed)//fathom on bounds
      milpBound = 1e200;
    }
#ifdef OA_DEBUG
    printf("Obj value after cuts %g %d rows\n",lp->getObjValue(),
	   numberCuts) ;
#endif
        //do we perform a new local search ?
        if (milpOptimal && feasible && !isInteger &&
            nLocalSearch_ < parameters_.maxLocalSearch_ &&
            numberPasses < parameters_.maxLocalSearchPerNode_ &&
            parameters_.localSearchNodeLimit_ > 0 &&
            CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_) {

            /** do we have a subMip? if not create a new one. */
            if(subMip == NULL) subMip = new SubMipSolver(lp, parameters_.strategy()); 

            nLocalSearch_++;

            subMip->performLocalSearch(cutoff, parameters_.subMilpLogLevel_,
                                       parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime(),
                                       parameters_.localSearchNodeLimit_);

            milpBound = subMip->lowBound();

            if(subMip->optimal())
            handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
            else
            handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;


            colsol = subMip->getLastSolution();
            isInteger = colsol != 0;

            feasible =  (milpBound < cutoff);

            if(feasible && isInteger)
             {
              bool changed = nlpManip.isDifferentOnIntegers(colsol);//If integer solution is the same as nlp
                                                                   //solution problem is solved
              if (!changed) {
                feasible = 0;
                milpBound = 1e50;
              }
              milpFeasible = feasible;
            }
            if (subMip->optimal())
              milpOptimal = 1;
            else {
              milpOptimal = 0;
            }

          if (milpBound < cutoff)
            handler_->message(UPDATE_LB, messages_)<<milpBound<<CoinCpuTime() - timeBegin_<<CoinMessageEol;
          else {
            milpBound = 1e50;
            feasible = 0;
            handler_->message(OASUCCESS, messages_)<<CoinCpuTime() - timeBegin_ <<CoinMessageEol;
          }
        }/** endif localSearch*/
    }
#ifdef OA_DEBUG
  debug_.printEndOfProcedureDebugMessage(cs, foundSolution, milpBound, isInteger, feasible, std::cout);
#endif

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
    if(leaveSiUnchanged_)
      lpManip->restore();
    delete lpManip;
    nlpManip.restore();
  return;
}


}/* End namespace Bonmin. */
