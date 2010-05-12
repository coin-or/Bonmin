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

#include <climits>

namespace Bonmin
{
   static const char * txt_id = "FP for MINLP";
/// Constructor with basic setup
  MinlpFeasPump::MinlpFeasPump(BabSetupBase & b):
      OaDecompositionBase(b, true, false)
  {
    int ivalue;
    std::string bonmin="bonmin.";
    std::string prefix = (b.prefix() == bonmin) ? "" : b.prefix();
    prefix += "pump_for_minlp.";
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
    parameter().maxLocalSearch_ = INT_MAX;
    b.options()->GetIntegerValue("solution_limit", parameter().maxSols_,prefix);
    parameter().maxLocalSearchTime_ =
    std::min(b.getDoubleParameter(BabSetupBase::MaxTime), oaTime);
    if(parameter().maxSols_ > b.getIntParameter(BabSetupBase::MaxSolutions))
      parameter().maxSols_ = b.getIntParameter(BabSetupBase::MaxSolutions);
       
  }

  MinlpFeasPump::~MinlpFeasPump()
  {}

  /// virutal method to decide if local search is performed
  bool
  MinlpFeasPump::doLocalSearch(BabInfo * babInfo) const
  {
    return (nLocalSearch_<parameters_.maxLocalSearch_ &&
        CoinCpuTime() - timeBegin_ < parameters_.maxLocalSearchTime_ &&
        numSols_ < parameters_.maxSols_);
  }
  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  double
  MinlpFeasPump::performOa(OsiCuts &cs,
      solverManip &lpManip,
      SubMipSolver * &subMip,
      BabInfo * babInfo,
      double & cutoff,const CglTreeInfo &info) const
  {

    //bool interuptOnLimit = false;
    //double lastPeriodicLog = CoinCpuTime();

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
    const double * colsol = NULL;
    lp->resolve();
    OsiBranchingInformation branch_info(lp, false);
    branch_info.lower_ = savedColLower();
    branch_info.upper_ = savedColUpper();
    if(lp->getNumCols() == nlp_->getNumCols())
      nlp_->addObjectiveFunction(*lp, nlp_->getColSolution());
    lp->setObjCoeff(numcols,0);

    bool milpOptimal = false;
    nlp_->resolve(txt_id);
    if (subMip)//Perform a local search
    {
      assert(subMip->solver() == lp);
      set_fp_objective(*lp, nlp_->getColSolution());
      lp->initialSolve();
      lp->setColUpper(numcols, cutoff);
      subMip->find_good_sol(DBL_MAX, parameters_.subMilpLogLevel_,
      //subMip->optimize(DBL_MAX, parameters_.subMilpLogLevel_,
          (parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime()) /* time limit */);

      milpOptimal = subMip -> optimal(); 
      colsol = subMip->getLastSolution();
      nLocalSearch_++;
      if(milpOptimal)
        handler_->message(SOLVED_LOCAL_SEARCH, messages_)
        <<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
      else
        handler_->message(LOCAL_SEARCH_ABORT, messages_)
        <<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
    }
    else
      colsol = lp->getColSolution();
    int numberPasses = 0;

#ifdef OA_DEBUG
    bool foundSolution = 0;
#endif
    double * nlpSol = NULL;
    int major_iteration = 0;
    double ub = cutoff;
    while (colsol) {
      numberPasses++;

      //after a prescribed elapsed time give some branch_information to user
      //double time = CoinCpuTime();


      //setup the nlp
      int numberCutsBefore = cs.sizeRowCuts();

      //Fix the variable which have to be fixed, after having saved the bounds
      branch_info.solution_ = colsol;

      vector<double> x_bar(indices.size());
      for(unsigned int i = 0 ; i < indices.size() ; i++){
         assert(fabs(colsol[indices[i]] - floor(colsol[indices[i]] + 0.5)) < 1e-5);
         x_bar[i] = colsol[indices[i]];
      }

      double dist = nlp_->solveFeasibilityProblem(indices.size(), x_bar(), indices(), 1, 0, 2);

      handler_->message(FP_DISTANCE, messages_) 
      <<dist<<CoinMessageEol;

      if(dist < 1e-05){
         fixIntegers(*nlp_,branch_info, parameters_.cbcIntegerTolerance_, objects_, nObjects_);

         nlp_->resolve(txt_id);
         if(!nlp_->isProvenOptimal()){
           relaxIntegers(*nlp_,branch_info, parameters_.cbcIntegerTolerance_, objects_, nObjects_);
           nlp_->resolve(txt_id);
         }
         bool restart = false;
         if (post_nlp_solve(babInfo, cutoff)) {
           restart = true;
           //nlp is solved and feasible
           // Update the cutoff
           ub = std::min(ub, nlp_->getObjValue());
           cutoff = ub * (1 - parameters_.cbcCutoffIncrement_);
           
           numSols_++;
         }
         else{
           //nlp_->setColLower(savedColLower());
           //nlp_->setColUpper(savedColUpper());
           //dist = nlp_->solveFeasibilityProblem(indices.size(), x_bar(), indices(), 1, 0, 2);
         }
         nlpSol = const_cast<double *>(nlp_->getColSolution());
         nlp_->getOuterApproximation(cs, nlpSol, 1, NULL,
                                  parameter().global_);
         //if(restart){
           nlp_->setColLower(savedColLower());
           nlp_->setColUpper(savedColUpper());
        if(restart){
          major_iteration++;
          handler_->message(FP_MAJOR_ITERATION, messages_) 
          <<major_iteration<<cutoff<<CoinMessageEol;
          nlp_->resolve(txt_id);
        }

         //}
      }
      else {
         nlpSol = const_cast<double *>(nlp_->getColSolution());
         nlp_->getOuterApproximation(cs, nlpSol, 1, NULL,
                                  parameter().global_);
      }


#if 0
      handler_->message(FP_MINOR_ITERATION, messages_) 
      <<nLocalSearch_<<cutoff<<CoinMessageEol;
#endif
      
      int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
      assert(numberCuts);
      installCuts(*lp, cs, numberCuts);
      numberCutsBefore = cs.sizeRowCuts();
     
      //check time
      if (CoinCpuTime() - timeBegin_ > parameters_.maxLocalSearchTime_)
        break;
      //do we perform a new local search ?
      if (nLocalSearch_ < parameters_.maxLocalSearch_ &&
          numSols_ < parameters_.maxSols_) {

        /** do we have a subMip? if not create a new one. */
        if (subMip == NULL) subMip = new SubMipSolver(lp, parameters_.strategy());

        nLocalSearch_++;
        set_fp_objective(*lp, nlp_->getColSolution());

        lp->setColUpper(numcols, cutoff); 
     
        subMip->find_good_sol(DBL_MAX, parameters_.subMilpLogLevel_,
        // subMip->optimize(DBL_MAX, parameters_.subMilpLogLevel_,
                         parameters_.maxLocalSearchTime_ + timeBegin_ - CoinCpuTime());
        milpOptimal = subMip -> optimal(); 
        colsol = subMip->getLastSolution();
      if(milpOptimal)
        handler_->message(SOLVED_LOCAL_SEARCH, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
      else
        handler_->message(LOCAL_SEARCH_ABORT, messages_)<<subMip->nodeCount()<<subMip->iterationCount()<<CoinMessageEol;
      if(colsol)
        handler_->message(FP_MILP_VAL, messages_) 
        <<colsol[nlp_->getNumCols()]<<CoinMessageEol;
         
      }/** endif localSearch*/
      else if (subMip!=NULL) {
        delete subMip;
        subMip = NULL;
        colsol = NULL;
        milpOptimal = false;
      }
    }
    if(colsol || ! milpOptimal)
      return -DBL_MAX;
    else{
      handler_->message(OASUCCESS, messages_)<<"FP"<<CoinCpuTime() - timeBegin_ 
      <<ub<<CoinMessageEol;
      return DBL_MAX;
    }
  }

  /** Register OA options.*/
  void
  MinlpFeasPump::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Options for feasibility pump", RegisteredOptions::BonminCategory);

    roptions->AddBoundedIntegerOption("fp_log_level",
        "specify FP iterations log level.",
        0,2,1,
        "Set the level of output of OA decomposition solver : "
        "0 - none, 1 - normal, 2 - verbose"
                                     );
    roptions->setOptionExtraInfo("fp_log_level",3);

    roptions->AddLowerBoundedNumberOption("fp_log_frequency",
        "display an update on lower and upper bounds in FP every n seconds",
        0.,1.,100.,
        "");
    roptions->setOptionExtraInfo("fp_log_frequency",3);
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
