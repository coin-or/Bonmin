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
    vector<double> savedColLower(nlp_->getNumCols());
    CoinCopyN(nlp_->getColLower(), nlp_->getNumCols(), savedColLower());
    vector<double> savedColUpper(nlp_->getNumCols());
    CoinCopyN(nlp_->getColUpper(), nlp_->getNumCols(), savedColUpper());


    OsiSolverInterface * lp = lpManip.si();

    vector<int> indices;
    for(int i = 0; i < numcols ; i++) {
      lp->setObjCoeff(i,0);
      if(!lp->isInteger(i)) {
      }
      else { indices.push_back(i);}
    }

    // If objective is linear need to add to lp constraint for objective
    if(lp->getNumCols() == nlp_->getNumCols())
      nlp_->addObjectiveFunction(*lp, nlp_->getColSolution());
    lp->setObjCoeff(numcols,0);
    const double * colsol = NULL;
    OsiBranchingInformation info(lp, false);

    nlp_->resolve();
    if (subMip)//Perform a local search
    {
      assert(subMip->solver() == lp);
      set_fp_objective(*lp, nlp_->getColSolution());
      lp->initialSolve();
      printf("Objective value %g\n", lp->getObjValue());
      lp->setColUpper(numcols, cutoff);
      subMip->find_good_sol(DBL_MAX, parameters_.subMilpLogLevel_,
          (parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime()) /* time limit */,
          parameters_.localSearchNodeLimit_);

      colsol = subMip->getLastSolution();
      nLocalSearch_++;

    }
    int numberPasses = 0;

#ifdef OA_DEBUG
    bool foundSolution = 0;
#endif
    double * nlpSol = NULL;

    while (colsol) {
      numberPasses++;

      //after a prescribed elapsed time give some information to user
      double time = CoinCpuTime();


      //setup the nlp
      int numberCutsBefore = cs.sizeRowCuts();

      //Fix the variable which have to be fixed, after having saved the bounds
      info.solution_ = colsol;

      vector<double> x_bar(indices.size());
      for(unsigned int i = 0 ; i < indices.size() ; i++){
         x_bar[i] = colsol[indices[i]];
      }

      double dist = nlp_->solveFeasibilityProblem(indices.size(), x_bar(), indices(), 1, 0, 2);

      printf("NLP solution is %g from MILP sol\n",dist);

      if(dist < 1e-06){
         fixIntegers(*nlp_,info, parameters_.cbcIntegerTolerance_, objects_, nObjects_);

         nlp_->resolve();
         if (post_nlp_solve(babInfo, cutoff)) {
           //nlp is solved and feasible
           // Update the cutoff
           cutoff = nlp_->getObjValue() * (1 - parameters_.cbcCutoffIncrement_);
           // Update the lp solver cutoff
           lp->setDblParam(OsiDualObjectiveLimit, cutoff);
           numSols_++;
         }
         nlpSol = const_cast<double *>(nlp_->getColSolution());
         nlp_->getOuterApproximation(cs, nlpSol, 1, NULL,
                                  parameter().global_);
         nlp_->setColLower(savedColLower());
         nlp_->setColUpper(savedColUpper());
         nlp_->resolve();
      }
      else {
         nlpSol = const_cast<double *>(nlp_->getColSolution());
         nlp_->getOuterApproximation(cs, nlpSol, 1, NULL,
                                  parameter().global_);
      }



      int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
      assert(numberCuts);
      installCuts(*lp, cs, numberCuts);
      numberCutsBefore = cs.sizeRowCuts();
     
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
        set_fp_objective(*lp, nlp_->getColSolution());

        lp->setColUpper(numcols, cutoff);

     
        subMip->find_good_sol(DBL_MAX, parameters_.subMilpLogLevel_,
            parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime(),
            parameters_.localSearchNodeLimit_);

        colsol = subMip->getLastSolution();
      }/** endif localSearch*/
      else if (subMip!=NULL) {
        delete subMip;
        subMip = NULL;
        colsol = NULL;
      }
    }
    if(colsol)
      return -DBL_MAX;
    else
      return DBL_MAX;
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
  si.initialSolve();
}

}/* End namespace Bonmin. */
