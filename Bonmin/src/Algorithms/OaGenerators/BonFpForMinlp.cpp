// (C) Copyright CNRS 2008
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, CNRS
//
// Date : 02/13/2009

#include "BonFpForMinlp.hpp"
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
  MinlpFeasPump::MinlpFeasPump(BabSetupBase & b):
      OaDecompositionBase(b, true, false)
  {
    int ivalue;
    b.options()->GetEnumValue("milp_subsolver",ivalue,b.prefix());
    if (ivalue <= 0) {//uses cbc
      //nothing to do?
    }
    else if (ivalue == 1) {
      int nodeS, nStrong, nTrust, mig, mir, probe, cover;
      b.options()->GetEnumValue("node_comparison",nodeS,"fp_sub.");
      b.options()->GetIntegerValue("number_strong_branch",nStrong,"fp_sub.");
      b.options()->GetIntegerValue("number_before_trust", nTrust,"fp_sub.");
      b.options()->GetIntegerValue("Gomory_cuts", mig,"fp_sub.");
      b.options()->GetIntegerValue("probing_cuts",probe,"fp_sub.");
      b.options()->GetIntegerValue("mir_cuts",mir,"fp_sub.");
      b.options()->GetIntegerValue("cover_cuts",cover,"fp_sub.");
      
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
    b.options()->GetNumericValue("minlp_pump_time_limit",oaTime,b.prefix());
    parameter().localSearchNodeLimit_ = 1000000;
    parameter().maxLocalSearch_ = 100000;
    parameter().maxLocalSearchPerNode_ = 10000;
    parameter().maxLocalSearchTime_ =
      Ipopt::Min(b.getDoubleParameter(BabSetupBase::MaxTime), oaTime);
  }
  MinlpFeasPump::~MinlpFeasPump()
  {}

  /// virutal method to decide if local search is performed
  bool
  MinlpFeasPump::doLocalSearch(BabInfo * babInfo) const
  {
    return (nLocalSearch_<parameters_.maxLocalSearch_ &&
        parameters_.localSearchNodeLimit_ > 0 &&
        CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_ &&
        numSols_ < parameters_.maxSols_);
  }
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  double
  MinlpFeasPump::performOa(OsiCuts &cs,
      solverManip &lpManip,
      SubMipSolver * &subMip,
      BabInfo * babInfo,
      double & cutoff) const
  {

    bool interuptOnLimit = false;
    double lastPeriodicLog= CoinCpuTime();

    const int numcols = nlp_->getNumCols();


    bool isInteger = false;

    OsiSolverInterface * lp = lpManip.si();
    OsiBranchingInformation info(lp, false);
    bool milpOptimal = 1;


    double milpBound = -COIN_DBL_MAX;
    bool milpFeasible = 1;
    bool feasible = 1;

    vector<int> indices;
    for(int i = 0; i < numcols ; i++) {
      if(!lp->isInteger(i)) {
         lp->setObjCoeff(i,0);
      }
      else { indices.push_back(i);}
    }

    if (subMip)//Perform a local search
    {
      set_fp_objective(*lp, nlp_->getColSolution());
      subMip->find_good_sol(cutoff, parameters_.subMilpLogLevel_,
          (parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime()) /* time limit */,
          parameters_.localSearchNodeLimit_);

      feasible = milpBound < cutoff;
      milpFeasible = feasible;
      isInteger = subMip->getLastSolution() != NULL;
      nLocalSearch_++;

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


      //setup the nlp
      int numberCutsBefore = cs.sizeRowCuts();

      //Fix the variable which have to be fixed, after having saved the bounds
      const double * colsol = subMip == NULL ? lp->getColSolution():
          subMip->getLastSolution();
      info.solution_ = colsol;

      vector<double> x_bar(indices.size());
      for(unsigned int i = 0 ; i < indices.size() ; i++){
         x_bar[i] = colsol[indices[i]];
      }
      nlp_->solveFeasibilityProblem(indices.size(), x_bar(), indices(), 1, 0, 2);

      info.solution_ = nlp_->getColSolution();

      if(integerFeasible(*lp, info, parameters_.cbcIntegerTolerance_,
                         objects_, nObjects_)){
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
       break;
      }


      nlpSol = const_cast<double *>(nlp_->getColSolution());

      // Get the cuts outer approximation at the current point
      const double * toCut = (parameter().addOnlyViolated_)?
          colsol:NULL;
      nlp_->getOuterApproximation(cs, nlpSol, 1, toCut,
                                  parameter().global_);

      int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
      if (numberCuts > 0)
        installCuts(*lp, cs, numberCuts);

      //check time
      if (CoinCpuTime() - timeBegin_ > parameters_.maxLocalSearchTime_)
        break;
      //do we perform a new local search ?
      if (nLocalSearch_ < parameters_.maxLocalSearch_ &&
          numberPasses < parameters_.maxLocalSearchPerNode_ &&
          parameters_.localSearchNodeLimit_ > 0 && numSols_ < parameters_.maxSols_) {

        /** do we have a subMip? if not create a new one. */
        if (subMip == NULL) subMip = new SubMipSolver(lp, parameters_.strategy());

        nLocalSearch_++;
        set_fp_objective(*lp_, nlp_->getColSolution());

     
        subMip->find_good_sol(cutoff, parameters_.subMilpLogLevel_,
            parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime(),
            parameters_.localSearchNodeLimit_);

        milpBound = subMip->lowBound();

        colsol = const_cast<double *> (subMip->getLastSolution());

      }/** endif localSearch*/
      else if (subMip!=NULL) {
        delete subMip;
        subMip = NULL;
      }
    }
    return -DBL_MAX;
  }

  /** Register OA options.*/
  void
  MinlpFeasPump::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Options for feasibility pump", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedNumberOption("minlp_pump_time_limit",
        "Specify the maximum number of seconds spent overall in MINLP Feasibility Pump.",
        0.,0,60.,
        "");

    roptions->AddBoundedIntegerOption("fp_log_level",
        "specify OA iterations log level.",
        0,2,1,
        "Set the level of output of OA decomposition solver : "
        "0 - none, 1 - normal, 2 - verbose"
                                     );

    roptions->AddLowerBoundedNumberOption("fp_log_frequency",
        "display an update on lower and upper bounds in OA every n seconds",
        0.,1.,100.,
        "");
  }

/** Put objective of MIP according to FP scheme. */
void
MinlpFeasPump::set_fp_objective(OsiSolverInterface &si, const double * colsol) const{
  if (objects_) {
    for (int i = 0 ; i < nObjects_ ; i++) {
      OsiObject * obj = objects_[i];
      int colnum = obj->columnNumber();
      if (colnum >= 0) {//Variable branching object
        double round = floor(colsol[colnum] + 0.5);
        double coeff = (colsol[colnum] - round ) < 0;
        si.setObjCoeff(colnum, 1 - 2 * coeff);
      }
      else {
        throw CoinError("OaDecompositionBase::solverManip",
                        "setFpObjective",
                        "Can not use FP on problem with SOS constraints");
      }
    }
  }
  else {
    int numcols = nlp_->getNumCols();
    for (int i = 0; i < numcols ; i++) {
      if (nlp_->isInteger(i)){
         double round = floor(colsol[i] + 0.5);
         double coeff = (colsol[i] - round ) < 0;
         si.setObjCoeff(i, 1 - 2*coeff);
      }
    }
  }
}

}/* End namespace Bonmin. */
