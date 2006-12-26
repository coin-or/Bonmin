// (C) Copyright International Business Machines 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 12/26/2006
//#define OA_DEBUG

#include "BonOaFeasChecker.hpp"
#include "BonminConfig.h"

#include "OsiAuxInfo.hpp"



namespace Bonmin
{


/// Default constructor
  OaFeasibilityChecker ::OaFeasibilityChecker ():
      CglCutGenerator(),
      OaDecompositionHelper()
   {
   }



  OaFeasibilityChecker ::OaFeasibilityChecker 
  (OsiTMINLPInterface * nlp,
   OsiSolverInterface * si,
   double cbcCutoffIncrement,
   double cbcIntegerTolerance,
   bool leaveSiUnchanged
   )
      :
      CglCutGenerator(),
      OaDecompositionHelper(nlp,si,
                          NULL, cbcCutoffIncrement,
                          cbcIntegerTolerance, leaveSiUnchanged)
  {
  }

  OaFeasibilityChecker ::~OaFeasibilityChecker ()
  {
  }

  /// cut generation method
  void
  OaFeasibilityChecker ::generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
      const CglTreeInfo info) const
  {
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
        return;
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
          "OaFeasibilityChecker ","generateCuts");
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

    bool feasible = 1;

    
    double milpBound = cutoff;
    int numberPasses = 0;
    while (isInteger && feasible ) {
      numberPasses++;

    //setup the nlp
    int numberCutsBefore = cs.sizeRowCuts();

    //Fix the variable which have to be fixed, after having saved the bounds
    colsol = (subMip == NULL) ? lp->getColSolution() : subMip->getLastSolution();
    nlpManip.fixIntegers(colsol);


      //Now solve the NLP get the cuts, and intall them in the local LP

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

      //  Reset the two solvers
      if(leaveSiUnchanged_)
        lpManip->restore();
      delete lpManip;
      nlpManip.restore();
    return;
  }

}/* End namespace Bonmin. */
