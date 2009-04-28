// (C) Copyright Carnegie Mellon University 2005, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
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


namespace Bonmin
{

/// Constructor with basic setup
  OACutGenerator2::OACutGenerator2(BabSetupBase & b):
      OaDecompositionBase(b, true, false)
  {
    int ivalue;
    std::string bonmin="bonmin.";
    std::string prefix = (b.prefix() == bonmin) ? "" : b.prefix();
    prefix += "oa_decomposition.";
    b.options()->GetEnumValue("milp_solver",ivalue,prefix);
    if (ivalue <= 0) {//uses cbc
      CbcStrategyDefault strategy;
      setStrategy(strategy);
    }
    else if (ivalue == 1) {
      CbcStrategyChooseCuts strategy(b, prefix);
      setStrategy(strategy);
    }
    else if (ivalue == 2) {
#ifdef COIN_HAS_CPX
      OsiCpxSolverInterface * cpxSolver = new OsiCpxSolverInterface;
      b.nonlinearSolver()->extractLinearRelaxation(*cpxSolver);
      assignLpInterface(cpxSolver);
#else
      std::cerr	<< "You have set an option to use CPLEX as the milp\n"
      << "subsolver in oa decomposition. However, apparently\n"
      << "CPLEX is not configured to be used in bonmin.\n"
      << "See the manual for configuring CPLEX\n";
      throw -1;
#endif
    }

    double oaTime;
    b.options()->GetNumericValue("time_limit",oaTime,prefix);
    parameter().maxLocalSearchTime_ =
    std::min(b.getDoubleParameter(BabSetupBase::MaxTime), oaTime);
    parameter().localSearchNodeLimit_ = 1000000;
    parameter().maxLocalSearch_ = 100000;
    parameter().maxLocalSearchPerNode_ = 10000;
  }
  OACutGenerator2::~OACutGenerator2()
  {}

  /// virutal method to decide if local search is performed
  bool
  OACutGenerator2::doLocalSearch(BabInfo * babInfo) const
  {
    return (nLocalSearch_<parameters_.maxLocalSearch_ &&
        parameters_.localSearchNodeLimit_ > 0 &&
        CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_ &&
        numSols_ < parameters_.maxSols_);
  }
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  double
  OACutGenerator2::performOa(OsiCuts &cs,
      solverManip &lpManip,
      SubMipSolver * &subMip,
      BabInfo * babInfo,
      double & cutoff) const
  {

    bool interuptOnLimit = false;
    double lastPeriodicLog = CoinCpuTime();

    //const int numcols = nlp_->getNumCols();


    bool isInteger = false;

    OsiSolverInterface * lp = lpManip.si();
    OsiBranchingInformation info(lp, false);
    bool milpOptimal = 1;


    double milpBound = -COIN_DBL_MAX;
    bool milpFeasible = 1;
    bool feasible = 1;

    if (subMip)//Perform a local search
    {
      subMip->find_good_sol(cutoff, parameters_.subMilpLogLevel_,
          (parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime()) /* time limit */,
          parameters_.localSearchNodeLimit_);
      milpBound = subMip->lowBound();
      milpOptimal = subMip->optimal();

      feasible = milpBound < cutoff;
      milpFeasible = feasible;
      isInteger = (subMip->getLastSolution() != NULL);
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
    double * nlpSol = NULL;

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
      const double * colsol = subMip == NULL ? lp->getColSolution():
          subMip->getLastSolution();
      info.solution_ = colsol;

      fixIntegers(*nlp_,info, parameters_.cbcIntegerTolerance_,objects_, nObjects_);

      nlp_->resolve();
      if (post_nlp_solve(babInfo, cutoff)) {
        //nlp solved and feasible
        // Update the cutoff
        cutoff = nlp_->getObjValue() *(1 - parameters_.cbcCutoffIncrement_);
        // Update the lp solver cutoff
        lp->setDblParam(OsiDualObjectiveLimit, cutoff);
        numSols_++;
      }

      nlpSol = const_cast<double *>(nlp_->getColSolution());

      // Get the cuts outer approximation at the current point
      const double * toCut = (parameter().addOnlyViolated_)?
          colsol:NULL;
      nlp_->getOuterApproximation(cs, nlpSol, 1, toCut,
                                  parameter().global_);

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
      info.solution_ = lp->getColSolution();
      if (!changed)
        changed = isDifferentOnIntegers(*nlp_, objects_, nObjects_,
                                        0.1,
                                        nlp_->getColSolution(), lp->getColSolution());
      if (changed) {

        isInteger = integerFeasible(*lp, info, parameters_.cbcIntegerTolerance_,
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
      //do we perform a new local search ?
      if (feasible && !isInteger &&
          nLocalSearch_ < parameters_.maxLocalSearch_ &&
          numberPasses < parameters_.maxLocalSearchPerNode_ &&
          parameters_.localSearchNodeLimit_ > 0 && numSols_ < parameters_.maxSols_) {

        /** do we have a subMip? if not create a new one. */
        if (subMip == NULL) subMip = new SubMipSolver(lp, parameters_.strategy());

        nLocalSearch_++;

        subMip->find_good_sol(cutoff, parameters_.subMilpLogLevel_,
            parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime(),
            parameters_.localSearchNodeLimit_);

        milpBound = subMip->lowBound();

        if (subMip->optimal())
          handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
        else
          handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;


        colsol = const_cast<double *> (subMip->getLastSolution());
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
      else if (subMip!=NULL) {
        delete subMip;
        subMip = NULL;
      }
    }

#ifdef OA_DEBUG
    debug_.printEndOfProcedureDebugMessage(cs, foundSolution, cutoff, milpBound, isInteger, feasible, std::cout);
#endif
    return milpBound;
  }

  /** Register OA options.*/
  void
  OACutGenerator2::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Options for OA decomposition", RegisteredOptions::BonminCategory);
    roptions->AddStringOption2("oa_decomposition", "If yes do initial OA decomposition",
                               "no",
                               "no","",
                               "yes","",
                               "");
    roptions->AddBoundedIntegerOption("oa_log_level",
        "specify OA iterations log level.",
        0,2,1,
        "Set the level of output of OA decomposition solver : "
        "0 - none, 1 - normal, 2 - verbose"
                                     );

    roptions->AddLowerBoundedNumberOption("oa_log_frequency",
        "display an update on lower and upper bounds in OA every n seconds",
        0.,1.,100.,
        "");
  }
}/* End namespace Bonmin. */

