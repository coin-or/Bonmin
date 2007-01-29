// (C) Copyright Carnegie Mellon University 2005, 2006
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
extern int usingCouenne;

/// Default constructor
  OACutGenerator2::OACutGenerator2():
      OaDecompositionBase()
   {
   }



  OACutGenerator2::OACutGenerator2
  (OsiTMINLPInterface * nlp,
   OsiSolverInterface * si,
   CbcStrategy * strategy,
   double cbcCutoffIncrement,
   double cbcIntegerTolerance,
   bool leaveSiUnchanged
   )
      :
      OaDecompositionBase(nlp,si,
                          strategy, cbcCutoffIncrement,
                          cbcIntegerTolerance, leaveSiUnchanged)
  {
  }

  OACutGenerator2::~OACutGenerator2()
  {
  }

  /// virutal method to decide if local search is performed
  bool
  OACutGenerator2::doLocalSearch() const{
    return (nLocalSearch_<parameters_.maxLocalSearch_ &&
          parameters_.localSearchNodeLimit_ > 0 &&
          CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_);
  }
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  double
  OACutGenerator2::performOa(OsiCuts &cs,
                             solverManip &nlpManip, 
                             solverManip &lpManip, 
                             SubMipSolver * &subMip,
                             OsiBabSolver * babInfo,
                             double & cutoff) const
  {
    double lastPeriodicLog= CoinCpuTime();

    const int numcols = nlp_->getNumCols();


    bool isInteger = false;

    OsiSolverInterface * lp = lpManip.si();
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

#ifdef OA_DEBUG
    bool foundSolution = 0;
#endif

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
   double * colsol = const_cast<double *>(subMip == NULL ? lp->getColSolution():
					  subMip->getLastSolution());
   if(usingCouenne){
     colsol = new double [numcols];
     CoinCopyN(lp->getColSolution(), numcols, colsol);}
   nlpManip.fixIntegers(colsol);


    if(solveNlp(babInfo, cutoff)){
      //nlp solved and feasible
      // Update the cutoff
      cutoff = nlp_->getObjValue() *(1 - parameters_.cbcCutoffIncrement_);
      // Update the lp solver cutoff
      lp->setDblParam(OsiDualObjectiveLimit, cutoff);
    }

      
    // Get the cuts outer approximation at the current point
      if(usingCouenne)
        nlpManip.restore();
    nlp_->getOuterApproximation(cs);


    int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
    if (numberCuts > 0)
      lpManip.installCuts(cs, numberCuts);

    lp->resolve();
    double objvalue = lp->getObjValue();
    //milpBound = max(milpBound, lp->getObjValue());
    feasible = (lp->isProvenOptimal() &&
		!lp->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
    //if value of integers are unchanged then we have to get out
    bool changed = !feasible;//if lp is infeasible we don't have to check anything
    if(!changed)
	  if(usingCouenne)
	    changed = nlpManip.isDifferentOnIntegers(lp->getColSolution(),colsol);
	  else
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
    //check time
    if(CoinCpuTime() - timeBegin_ > parameters_.maxLocalSearchTime_)
      break;
        //do we perform a new local search ?
    if (milpOptimal && feasible && !isInteger &&
	nLocalSearch_ < parameters_.maxLocalSearch_ &&
	numberPasses < parameters_.maxLocalSearchPerNode_ &&
	parameters_.localSearchNodeLimit_ > 0) {
      
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


            colsol =const_cast<double *> (subMip->getLastSolution());
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
    else if(subMip!=NULL) {
      delete subMip; 
      subMip = NULL;
    }
    if(usingCouenne)
      delete [] colsol;
    }
#ifdef OA_DEBUG
  debug_.printEndOfProcedureDebugMessage(cs, foundSolution, milpBound, isInteger, feasible, std::cout);
#endif
  return milpBound;
}


}/* End namespace Bonmin. */
