// (C) Copyright Carnegie Mellon University 2005, 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 05/26/2005

#include "BonOACutGenerator2.hpp"
#include "BonminConfig.h"

#include "OsiClpSolverInterface.hpp"

#include "CbcModel.hpp"
#include "BonCbcLpStrategy.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#include "OsiAuxInfo.hpp"
#include "BonSolverHelp.hpp"

#include <climits>

namespace Bonmin
{
   static const char * txt_id = "OA decomposition";


/// Constructor with basic setup
  OACutGenerator2::OACutGenerator2(BabSetupBase & b):
      OaDecompositionBase(b, true, false)
  {
    std::string bonmin="bonmin.";
    std::string prefix = (b.prefix() == bonmin) ? "" : b.prefix();
    prefix += "oa_decomposition.";
    subMip_ = new SubMipSolver(b, prefix);
    double oaTime;
    b.options()->GetNumericValue("time_limit",oaTime,prefix);
    parameter().maxLocalSearch_ = INT_MAX;
    b.options()->GetIntegerValue("solution_limit", parameter().maxSols_,prefix);
    parameter().maxLocalSearchTime_ =
    std::min(b.getDoubleParameter(BabSetupBase::MaxTime), oaTime);
    if(parameter().maxSols_ > b.getIntParameter(BabSetupBase::MaxSolutions))
      parameter().maxSols_ = b.getIntParameter(BabSetupBase::MaxSolutions);
  }
  OACutGenerator2::~OACutGenerator2()
  {
     delete subMip_;
  }

  /// virutal method to decide if local search is performed
  bool
  OACutGenerator2::doLocalSearch(BabInfo * babInfo) const
  {
    return (nLocalSearch_<parameters_.maxLocalSearch_ &&
            numSols_ < parameters_.maxSols_ &&
	    CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_);
  }
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  double
  OACutGenerator2::performOa(OsiCuts &cs,
      solverManip &lpManip,
      BabInfo * babInfo,
      double & cutoff, const CglTreeInfo & info) const
  {

    double lastPeriodicLog = CoinCpuTime();

    //const int numcols = nlp_->getNumCols();
    double gap_tol = this->parameter().gap_tol_;

    bool isInteger = false;

    subMip_->setLpSolver(lpManip.si());
    OsiSolverInterface * lp = subMip_->solver();
    lp->resolve();

    if(IsValid(nlp_->linearizer())){
      nlp_->linearizer()->get_refined_oa(cs);
      installCuts(*lp, cs, cs.sizeRowCuts());
    }
    lp->resolve();

    OsiBranchingInformation branch_info(lp, false);
    bool milpOptimal = 1;


    double milpBound = -COIN_DBL_MAX;
    //bool milpFeasible = 1;
    bool feasible = 1;

    subMip_->solve(cutoff, parameters_.subMilpLogLevel_,
                   parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());

    milpBound = std::max(milpBound, subMip_->lowBound());
    milpOptimal = subMip_->optimal();

    feasible = milpBound < cutoff;
    //milpFeasible = feasible;
    isInteger = (subMip_->getLastSolution() != NULL);
    nLocalSearch_++;

    if (milpOptimal)
      handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip_->nodeCount()<<subMip_->iterationCount()<<CoinMessageEol;
    else
    {
      handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip_->nodeCount()<<subMip_->iterationCount()<<CoinMessageEol;
    }

    int nodeCount = subMip_->nodeCount();

    int numberPasses = 0;

#ifdef OA_DEBUG
    bool foundSolution = 0;
#endif
    double * nlpSol = NULL;
    double ub = cutoff;
    double gap = 1;
    while (isInteger && feasible) {
      numberPasses++;
      //after a prescribed elapsed time give some information to user
      double time = CoinCpuTime();
      if (time - lastPeriodicLog > parameters_.logFrequency_) {
        handler_->message(PERIODIC_MSG,messages_)
        <<time - timeBegin_<<numberPasses<<cutoff
        <<milpBound
        <<CoinMessageEol;
        lastPeriodicLog = CoinCpuTime();
      }


      //setup the nlp
      int numberCutsBefore = cs.sizeRowCuts();

      //Fix the variable which have to be fixed, after having saved the bounds
      const double * colsol = 
          subMip_->getLastSolution();
      branch_info.solution_ = colsol;

      fixIntegers(*nlp_,branch_info, parameters_.cbcIntegerTolerance_,objects_, nObjects_);

      nlp_->resolve(txt_id);
      if (post_nlp_solve(babInfo, cutoff)) {
        //nlp solved and feasible
        // Update the cutoff
        ub = std::min(nlp_->getObjValue(), ub);
        cutoff = ub > 0 ? ub *(1 - parameters_.cbcCutoffIncrement_) : ub*(1 + parameters_.cbcCutoffIncrement_);
        assert(cutoff < ub);
        // Update the lp solver cutoff
        lp->setDblParam(OsiDualObjectiveLimit, cutoff);
        numSols_++;
      }

      nlpSol = const_cast<double *>(nlp_->getColSolution());

      // Get the cuts outer approximation at the current point
      if(parameter().addOnlyViolated_){
        nlp_->getOuterApproximation(cs, nlpSol, 1, colsol,
                                     parameter().global_);
      }
      else {
        nlp_->getOuterApproximation(cs, nlpSol, 1, NULL,
                                     parameter().global_);
      }

      int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
      assert(numberCuts);
      installCuts(*lp, cs, numberCuts);

      lp->resolve();

      double objvalue = lp->getObjValue();
      //milpBound = max(milpBound, lp->getObjValue());
      feasible = (lp->isProvenOptimal() &&
          !lp->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
      //if value of integers are unchanged then we have to get out
      bool changed = !feasible;//if lp is infeasible we don't have to check anything
      branch_info.solution_ = lp->getColSolution();
      if (!changed)
        changed = isDifferentOnIntegers(*nlp_, objects_, nObjects_,
                                        parameters_.cbcIntegerTolerance_,
                                        nlp_->getColSolution(), lp->getColSolution());
      if (changed) {

        isInteger = integerFeasible(*lp, branch_info, parameters_.cbcIntegerTolerance_,
                                     objects_, nObjects_);
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
      if (CoinCpuTime() - timeBegin_ > parameters_.maxLocalSearchTime_)
        break;
      
      if(ub!=0)
    	  gap = (ub - milpBound)/fabs(ub);
      else
		  gap = -milpBound/(1e-10);
      if (gap < gap_tol){
		 milpBound = 1e50;
		 feasible = 0;
	   }
      //do we perform a new local search ?
      if (feasible && 
          nLocalSearch_ < parameters_.maxLocalSearch_ &&
	  numSols_ < parameters_.maxSols_) {

	  if(ub!=0)
		  gap = (ub - milpBound)/fabs(ub);
	  else
		  gap = -milpBound/(1e-10);
	  if (gap < gap_tol){
		 milpBound = 1e50;
		 feasible = 0;
	   }
        nLocalSearch_++;

        subMip_->solve(cutoff, parameters_.subMilpLogLevel_,
            parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());

        milpBound = std::max(milpBound, subMip_->lowBound());

        if (subMip_->optimal())
          handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip_->nodeCount()<<subMip_->iterationCount()<<CoinMessageEol;
        else
          handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip_->nodeCount()<<subMip_->iterationCount()<<CoinMessageEol;

        nodeCount += subMip_->nodeCount();

        colsol = const_cast<double *> (subMip_->getLastSolution());
        isInteger = (colsol != 0);

        feasible =  (milpBound < cutoff);

        if (feasible && isInteger) {
          bool changed = false;
          changed = isDifferentOnIntegers(*nlp_, objects_, nObjects_,
                                          0.1,
                                          nlp_->getColSolution(), colsol);
          //solution problem is solved
          if (!changed) {
            feasible = 0;
            milpBound = 1e50;
          }
          //milpFeasible = feasible;
        }
        if (subMip_->optimal())
          milpOptimal = 1;
        else {
          milpOptimal = 0;
        }
        if (milpBound < cutoff)
		  handler_->message(UPDATE_LB, messages_)<<milpBound<<CoinCpuTime() - timeBegin_<<CoinMessageEol;
		else {
		  milpBound = 1e50;
		  feasible = 0;
		}
      }
      
    }

   if(milpBound >= cutoff){
      handler_->message(OASUCCESS, messages_)<<"OA "<<CoinCpuTime() - timeBegin_ 
                                             <<ub<<milpBound<<CoinMessageEol;
   }
   else {
      handler_->message(OAABORT, messages_)<<"OA "<<CoinCpuTime() - timeBegin_           
                                             <<ub<<milpBound<<CoinMessageEol;
   }
   handler_->message(OA_STATS, messages_)<<numberPasses<<nodeCount
                                         <<CoinMessageEol;
#ifdef OA_DEBUG
    debug_.printEndOfProcedureDebugMessage(cs, foundSolution, cutoff, milpBound, isInteger, feasible, std::cout);
#endif
    return milpBound;
  }

  /** Register OA options.*/
  void
  OACutGenerator2::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Outer Approximation Decomposition (B-OA)", RegisteredOptions::BonminCategory);
    roptions->AddStringOption2("oa_decomposition", "If yes do initial OA decomposition",
                               "no",
                               "no","",
                               "yes","",
                               "");
    roptions->setOptionExtraInfo("oa_decomposition",19);

    roptions->SetRegisteringCategory("Output and Loglevel", RegisteredOptions::BonminCategory);
    roptions->AddBoundedIntegerOption("oa_log_level",
        "specify OA iterations log level.",
        0,2,1,
        "Set the level of output of OA decomposition solver : "
        "0 - none, 1 - normal, 2 - verbose"
                                     );
    roptions->setOptionExtraInfo("oa_log_level", 25);

    roptions->AddLowerBoundedNumberOption("oa_log_frequency",
        "display an update on lower and upper bounds in OA every n seconds",
        0.,1.,100.,
        "");
    roptions->setOptionExtraInfo("oa_log_frequency", 25);
  }
}/* End namespace Bonmin. */

