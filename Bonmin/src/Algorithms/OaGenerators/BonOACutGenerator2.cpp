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
#include "BonCbcLpStrategy.hpp"
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#include "OsiAuxInfo.hpp"

#include <climits>

namespace Bonmin
{

/// Constructor with basic setup
  OACutGenerator2::OACutGenerator2(BabSetupBase & b):
      OaDecompositionBase(b, true, false)
  {
    int ivalue;
    b.options()->GetEnumValue("milp_subsolver",ivalue,"bonmin.");
    if (ivalue <= 0) {//uses cbc
      //nothing to do?
    }
    else if (ivalue == 1) {
      int nodeS, nStrong, nTrust, mig, mir, probe, cover;
      b.options()->GetEnumValue("node_comparison",nodeS,"milp_sub.");
      b.options()->GetIntegerValue("number_strong_branch",nStrong,"milp_sub.");
      b.options()->GetIntegerValue("number_before_trust", nTrust,"milp_sub.");
      b.options()->GetIntegerValue("Gomory_cuts", mig,"milp_sub.");
      //b.options()->GetIntegerValue("probing_cuts",probe,"milp_sub.");
      probe = 0;
      b.options()->GetIntegerValue("mir_cuts",mir,"milp_sub.");
      b.options()->GetIntegerValue("cover_cuts",cover,"milp_sub.");
      
      CbcStrategy * strategy =
        new CbcOaStrategy(mig, probe, mir, cover, nTrust,
            nStrong, nodeS, parameters_.cbcIntegerTolerance_, parameters_.subMilpLogLevel_);
      setStrategy(*strategy);
      delete strategy;

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
    b.options()->GetNumericValue("oa_dec_time_limit",oaTime,"bonmin.");
    parameter().localSearchNodeLimit_ = INT_MAX;
    parameter().maxLocalSearch_ = INT_MAX;
    parameter().maxLocalSearchPerNode_ = INT_MAX;
    parameter().maxLocalSearchTime_ =
      Ipopt::Min(b.getDoubleParameter(BabSetupBase::MaxTime), oaTime);
  }
  OACutGenerator2::~OACutGenerator2()
  {}

  /// virutal method to decide if local search is performed
  bool
  OACutGenerator2::doLocalSearch() const
  {
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
      nlpManip.fixIntegers(info);


      if (solveNlp(babInfo, cutoff)) {
        //nlp solved and feasible
        // Update the cutoff
        cutoff = nlp_->getObjValue() *(1 - parameters_.cbcCutoffIncrement_);
        // Update the lp solver cutoff
        lp->setDblParam(OsiDualObjectiveLimit, cutoff);
      }

      nlpSol = const_cast<double *>(nlp_->getColSolution());

      // Get the cuts outer approximation at the current point
      const double * toCut = (parameter().addOnlyViolated_)?
          colsol:NULL;
      nlp_->getOuterApproximation(cs, nlpSol, 1, toCut,
                                  parameter().global_);

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
      info.solution_ = lp->getColSolution();
      if (!changed)
        changed = nlpManip.isDifferentOnIntegers(lp->getColSolution());
      if (changed) {

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
      //check time
      if (CoinCpuTime() - timeBegin_ > parameters_.maxLocalSearchTime_)
        break;
      //do we perform a new local search ?
      if (milpOptimal && feasible && !isInteger &&
          nLocalSearch_ < parameters_.maxLocalSearch_ &&
          numberPasses < parameters_.maxLocalSearchPerNode_ &&
          parameters_.localSearchNodeLimit_ > 0) {

        /** do we have a subMip? if not create a new one. */
        if (subMip == NULL) subMip = new SubMipSolver(lp, parameters_.strategy());

        nLocalSearch_++;

        subMip->performLocalSearch(cutoff, parameters_.subMilpLogLevel_,
            parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime(),
            parameters_.localSearchNodeLimit_);

        milpBound = subMip->lowBound();

        if (subMip->optimal())
          handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
        else
          handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;


        colsol = const_cast<double *> (subMip->getLastSolution());
        isInteger = colsol != 0;

        feasible =  (milpBound < cutoff);

        if (feasible && isInteger) {
          bool changed = false;
          changed = nlpManip.isDifferentOnIntegers(colsol);//If integer solution is the same as nlp
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
    roptions->AddLowerBoundedNumberOption("oa_dec_time_limit",
        "Specify the maximum number of seconds spent overall in OA decomposition iterations.",
        0.,0,30.,
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

    roptions->SetRegisteringCategory("Options for MILP subsolver in OA decomposition", RegisteredOptions::BonminCategory);
    roptions->AddStringOption3("milp_subsolver",
        "Choose the subsolver to solve MILP sub-problems in OA decompositions.",
        "Cbc_D",
        "Cbc_D","Coin Branch and Cut with its default",
        "Cbc_Par", "Coin Branch and Cut with passed parameters",
        "Cplex","Ilog Cplex",
        " To use Cplex, a valid license is required and you should have compiled OsiCpx in COIN-OR  (see Osi documentation).");
    roptions->setOptionExtraInfo("milp_subsolver",5);

    roptions->AddBoundedIntegerOption("milp_log_level",
        "specify MILP subsolver log level.",
        0,3,0,
        "Set the level of output of the MILP subsolver in OA : "
        "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                     );
    roptions->setOptionExtraInfo("milp_log_level",5);


  }


}/* End namespace Bonmin. */
