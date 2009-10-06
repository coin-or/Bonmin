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
//#define OA_DEBUG


namespace Bonmin
{

// Default constructor
  OaFeasibilityChecker ::OaFeasibilityChecker
  (OsiTMINLPInterface * nlp,
   OsiSolverInterface * si,
   double cbcCutoffIncrement,
   double cbcIntegerTolerance,
   bool leaveSiUnchanged
  )
      :
      OaDecompositionBase(nlp,si,
          NULL, cbcCutoffIncrement,
          cbcIntegerTolerance, leaveSiUnchanged)
  {}

  /// New usefull constructor
  OaFeasibilityChecker::OaFeasibilityChecker(BabSetupBase &b):
      OaDecompositionBase(b, false, true)
  {}
  OaFeasibilityChecker ::~OaFeasibilityChecker ()
  {}

  /// OaDecomposition method
  double
  OaFeasibilityChecker::performOa(OsiCuts & cs, solverManip &nlpManip, solverManip &lpManip,
      SubMipSolver * &subMip, OsiBabSolver * babInfo, double &cutoff) const
  {
    bool isInteger = true;
    bool feasible = 1;

    OsiCuts cs2;
    OsiSolverInterface * lp = lpManip.si();
    OsiBranchingInformation info(lp,false);
    //int numcols = lp->getNumCols();
    double milpBound = -COIN_DBL_MAX;
    int numberPasses = 0;
    double * nlpSol =  NULL;
    while (isInteger && feasible ) {
      numberPasses++;

      //setup the nlp
      int numberCutsBefore = cs2.sizeRowCuts();

      //Fix the variable which have to be fixed, after having saved the bounds
      double * colsol = const_cast<double *>(lp->getColSolution());
      info.solution_ = colsol;
      nlpManip.fixIntegers(info);


      //Now solve the NLP get the cuts, and intall them in the local LP

      if (solveNlp(babInfo, cutoff)) {
        //nlp solved and feasible
        // Update the cutoff
        cutoff = nlp_->getObjValue() *(1 - parameters_.cbcCutoffIncrement_);
        // Update the lp solver cutoff
        lp->setDblParam(OsiDualObjectiveLimit, cutoff);
      }
      // Get the cuts outer approximation at the current point

      nlpSol = const_cast<double *>(nlp_->getColSolution());

      const double * toCut = (parameter().addOnlyViolated_)?
          colsol:NULL;
      nlp_->getOuterApproximation(cs2, nlpSol, 1, toCut,
          parameter().global_);
      int numberCuts = cs2.sizeRowCuts() - numberCutsBefore;
      if (numberCuts > 0)
        lpManip.installCuts(cs2, numberCuts);

      lp->resolve();
      double objvalue = lp->getObjValue();
      //milpBound = max(milpBound, lp->getObjValue());
      feasible = (lp->isProvenOptimal() &&
          !lp->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
      //if value of integers are unchanged then we have to get out
      bool changed = true;//if lp is infeasible we don't have to check anything
      isInteger = 0;
      //	  if(!fixed)//fathom on bounds
      //           milpBound = 1e200;
      if (feasible) {
        changed = nlpManip.isDifferentOnIntegers(lp->getColSolution());
      }
      if (changed) {
        info.solution_ = lp->getColSolution();
        isInteger = lpManip.integerFeasible(info);
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
    int numberCuts = cs2.sizeRowCuts();
    //Remove non tight cuts
    if(numberCuts)
    {
       int num_rows_before = lp->getNumRows() - numberCuts;
       vector<int> toDelete;
       toDelete.reserve(numberCuts);
       CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (lp->getWarmStart());
       assert(basis);
       for(int i = 0 ; i < numberCuts ; i++){
          int idx = i + num_rows_before;
          if(basis->getArtifStatus(idx) == CoinWarmStartBasis::basic){
             toDelete.push_back(idx);
          }
          else {
            cs.insert(cs2.rowCut(i));
          }
       }
       lp->deleteRows(toDelete.size(), toDelete());
       //printf("Remove %i cuts\n", toDelete.size());
       lp->resolve();
       delete basis;
    } 
#ifdef OA_DEBUG
    //debug_.printEndOfProcedureDebugMessage(cs, foundSolution, milpBound, isInteger, feasible, std::cout);
    std::cout<<"milpBound found: "<<milpBound<<std::endl;
#endif
    return milpBound;
  }

}/* End namespace Bonmin. */
