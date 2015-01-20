// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006


// Code separated from BonOaDecBase to try to clarify OAs
#include "BonSubMipSolver.hpp"
#include "BonminConfig.h"
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiClpSolverInterface.hpp"

#include <climits>
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
void throw_error(const std::string &s, const std::string &f, const std::string &func){
throw CoinError(s,f,func);
}
#define CHECK_CPX_STAT(a,b) if(b) throw_error("Error in CPLEX call",__FILE__,a);

#endif

#include "BonRegisteredOptions.hpp"
#include "BonBabSetupBase.hpp"
#include "BonCbcLpStrategy.hpp"


namespace Bonmin {
  /** Constructor */
  SubMipSolver::SubMipSolver(BabSetupBase &b, const std::string &prefix):
      clp_(NULL),
      cpx_(NULL),
      lowBound_(-DBL_MAX),
      optimal_(false),
      integerSolution_(NULL),
      strategy_(NULL),
      ownClp_(false)
  {

   int logLevel;
   b.options()->GetIntegerValue("milp_log_level", logLevel, prefix);

   int ivalue;
   b.options()->GetEnumValue("milp_solver",ivalue,prefix);
   if (ivalue <= 0) {//uses cbc
     strategy_ = new CbcStrategyDefault;
     clp_ = new OsiClpSolverInterface;
     ownClp_ = true;
     clp_->messageHandler()->setLogLevel(logLevel);
   }
   else if (ivalue == 1) {
     CbcStrategyChooseCuts strategy(b, prefix);
     strategy_  = new CbcStrategyChooseCuts(b, prefix);
     clp_ = new OsiClpSolverInterface;
     ownClp_ = true;
     clp_->messageHandler()->setLogLevel(logLevel);
   }
   else if (ivalue == 2) {
#ifdef COIN_HAS_CPX
      OsiCpxSolverInterface * cpxSolver = new OsiCpxSolverInterface;
#if 1
      
      b.options()->GetIntegerValue("number_cpx_threads",ivalue,prefix);
      CPXsetintparam(cpxSolver->getEnvironmentPtr(), CPX_PARAM_THREADS, ivalue);
      b.options()->GetIntegerValue("cpx_parallel_strategy",ivalue,prefix);
      CPXsetintparam(cpxSolver->getEnvironmentPtr(), CPX_PARAM_PARALLELMODE, ivalue);
#endif
      cpx_ = cpxSolver;
      cpx_->messageHandler()->setLogLevel(logLevel);
#else
      std::cerr	<< "You have set an option to use CPLEX as the milp\n"
      << "subsolver in oa decomposition. However, apparently\n"
      << "CPLEX is not configured to be used in bonmin.\n"
      << "See the manual for configuring CPLEX\n";
      throw -1;
#endif
    }

      b.options()->GetEnumValue("milp_strategy",ivalue,prefix);
      if(ivalue == 0){
        milp_strat_ = FindGoodSolution;
      }
      else {
        milp_strat_ = GetOptimum;
      }

      b.options()->GetNumericValue("allowable_fraction_gap", gap_tol_, prefix);


  }
  SubMipSolver::SubMipSolver(const SubMipSolver &copy):
      clp_(NULL),
      cpx_(NULL),
      lowBound_(-DBL_MAX),
      optimal_(false),
      integerSolution_(NULL),
      strategy_(NULL),
      milp_strat_(copy.milp_strat_),
      gap_tol_(copy.gap_tol_),
      ownClp_(copy.ownClp_)
  {
#ifdef COIN_HAS_CPX
     if(copy.cpx_ != NULL){
       cpx_ = new OsiCpxSolverInterface(*copy.cpx_);
      int ival;
      CPXgetintparam(copy.cpx_->getEnvironmentPtr(), CPX_PARAM_THREADS, &ival);
      CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_THREADS, ival);
      CPXgetintparam(copy.cpx_->getEnvironmentPtr(), CPX_PARAM_PARALLELMODE, &ival);
      CPXsetintparam(cpx_->getEnvironmentPtr(), CPX_PARAM_PARALLELMODE, ival);
     }
#endif
     if(copy.clp_ != NULL){
       if(ownClp_) clp_ = new OsiClpSolverInterface(*copy.clp_);
       else clp_ = copy.clp_;
     }
     if(copy.strategy_){
        strategy_ = dynamic_cast<CbcStrategyDefault *>(copy.strategy_->clone());
        assert(strategy_);
     }
  }
  SubMipSolver::~SubMipSolver()
  {
    if (strategy_) delete strategy_;
    if (integerSolution_) delete [] integerSolution_;
    #ifdef COIN_HAS_CPX
    if(cpx_) delete cpx_;
    #endif
    if(ownClp_) delete clp_;
  }

  /** Assign lp solver. */
  void
  SubMipSolver::setLpSolver(OsiSolverInterface * lp)
  {
#ifdef COIN_HAS_CPX
    if(cpx_){
      clp_ = NULL;
      cpx_->loadProblem(*lp->getMatrixByCol(), lp->getColLower(), lp->getColUpper(), lp->getObjCoefficients(), lp->getRowLower(), lp->getRowUpper());
      int ncols = lp->getNumCols();
      for(int i = 0 ; i < ncols ; i++){
        if(lp->isInteger(i) || lp->isBinary(i))
           cpx_->setInteger(i);
        else
           cpx_->setContinuous(i);
      }
    }
    else {
#endif
      if(ownClp_) delete clp_;
      ownClp_ = false;
      clp_ = (lp == NULL) ? NULL :
              dynamic_cast<OsiClpSolverInterface *>(lp);
      assert(clp_);
#ifdef COIN_HAS_CPX
    }
#endif
    lowBound_ = -COIN_DBL_MAX;
    optimal_ = false;
    if (integerSolution_) {
      delete [] integerSolution_;
      integerSolution_ = NULL;
    }
  }

  OsiSolverInterface * 
  SubMipSolver::solver(){
         if(clp_ != NULL)
           return clp_;
         else
#ifdef COIN_HAS_CPX
           return cpx_;
#else
         return NULL;
#endif
      }

 void 
 SubMipSolver::find_good_sol(double cutoff, int loglevel, double max_time){

     if(clp_){
      CbcStrategyDefault * strat_default = NULL;
      if (!strategy_){
        strat_default = new CbcStrategyDefault(1,5,5, loglevel);
        strat_default->setupPreProcessing();
        strategy_ = strat_default;
      }
      OsiBabSolver empty;
      CbcModel cbc(*clp_);
      cbc.solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(cbc.messagesPointer()->source_,"OCbc");

      cbc.setLogLevel(loglevel);
      cbc.solver()->messageHandler()->setLogLevel(0);
      clp_->resolve();
      cbc.setStrategy(*strategy_);
      cbc.setLogLevel(loglevel);
      cbc.solver()->messageHandler()->setLogLevel(0);
      cbc.setMaximumSeconds(max_time);
      cbc.setMaximumSolutions(1);
      cbc.setCutoff(cutoff);

      
      cbc.branchAndBound();
      lowBound_ = cbc.getBestPossibleObjValue();

      if (cbc.isProvenOptimal() || cbc.isProvenInfeasible())
        optimal_ = true;
      else optimal_ = false;

      if (cbc.getSolutionCount()) {
        if (!integerSolution_)
          integerSolution_ = new double[clp_->getNumCols()];
        CoinCopyN(cbc.bestSolution(), clp_->getNumCols(), integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
      nodeCount_ = cbc.getNodeCount();
      iterationCount_ = cbc.getIterationCount();

      if(strat_default != NULL){
        delete strat_default;
        strategy_ = NULL;
      }
     }
     else if (cpx_){
#ifndef COIN_HAS_CPX
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
#else
        cpx_->messageHandler()->setLogLevel(loglevel);
        cpx_->switchToMIP();
        CPXENVptr env = cpx_->getEnvironmentPtr();
        CPXLPptr cpxlp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
        CPXsetdblparam(env, CPX_PARAM_TILIM, max_time);
        CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 1);
        CPXsetdblparam(env, CPX_PARAM_EPINT, 1e-08);
        CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);
        CPXsetdblparam(env, CPX_PARAM_EPGAP, gap_tol_);

        double start_time = CoinCpuTime();

        CPXsetintparam(env,CPX_PARAM_INTSOLLIM, 10);
        CPXsetintparam(env,CPX_PARAM_NODELIM, 1000);

        nodeCount_ = 0;
        iterationCount_ = 0;
        int status = CPXmipopt(env,cpxlp);
        CHECK_CPX_STAT("mipopt",status)
       
    
        status = CPXgetbestobjval(env, cpxlp, &lowBound_);
        CHECK_CPX_STAT("bestobjvalue",status)
     
        int stat = CPXgetstat( env, cpxlp);
        optimal_ = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL) || (stat == CPXMIP_INForUNBD) ; 
        nodeCount_ = CPXgetnodecnt(env , cpxlp);
        iterationCount_ = CPXgetmipitcnt(env , cpxlp);
        
        int type;
        status = CPXsolninfo(env, cpxlp, NULL, &type, NULL, NULL);
        CHECK_CPX_STAT("solninfo", status);
        while(!optimal_ && type == CPX_NO_SOLN && stat != CPXMIP_SOL_LIM && stat != CPXMIP_TIME_LIM_INFEAS 
              && stat != CPXMIP_TIME_LIM_FEAS && (CoinCpuTime() - start_time) <= max_time){
          CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
          CPXsetdblparam(env, CPX_PARAM_TILIM, max_time - CoinCpuTime() + start_time);
          CPXsetintparam(env,CPX_PARAM_NODELIM, 2100000000);
           status = CPXmipopt(env,cpxlp);
           CHECK_CPX_STAT("mipopt",status)

           stat = CPXgetstat( env, cpxlp);
           optimal_ = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL) || (stat == CPXMIP_INForUNBD) ; 
           nodeCount_ = CPXgetnodecnt(env , cpxlp);
           iterationCount_ = CPXgetmipitcnt(env , cpxlp);
       }
       bool infeasible = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_ABORT_INFEAS) || (stat == CPXMIP_TIME_LIM_INFEAS) 
                          || (stat == CPXMIP_NODE_LIM_INFEAS) || (stat == CPXMIP_FAIL_INFEAS)  
                          || (stat == CPXMIP_MEM_LIM_INFEAS) || (stat == CPXMIP_INForUNBD);
       

      status = CPXgetbestobjval(env, cpxlp, &lowBound_);
      CHECK_CPX_STAT("getbestobjval",status)
       if(!infeasible){
          nodeCount_ += CPXgetnodecnt(env, cpxlp);
          //iterationCount_ += CPXgetitcnt(env, cpxlp);
          if(!integerSolution_){
            integerSolution_ = new double[cpx_->getNumCols()];
          }
          CPXgetmipx(env, cpxlp, integerSolution_, 0, cpx_->getNumCols() -1);
          CHECK_CPX_STAT("getmipx",status)
       }
       else {
         if (integerSolution_) {
           delete [] integerSolution_;
           integerSolution_ = NULL;
         }
       }
      cpx_->switchToLP();
#endif
    } 
    else {
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
   }
  }

  void
  SubMipSolver::optimize(double cutoff, int loglevel, double maxTime)
  {
    if (clp_) {
      assert(strategy_);
      CbcStrategyDefault * strat_default = dynamic_cast<CbcStrategyDefault *>(strategy_->clone());
      assert(strat_default);
      strat_default->setupPreProcessing();

      OsiBabSolver empty;
      CbcModel cbc(*clp_);
      cbc.solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(cbc.messagesPointer()->source_,"OCbc");

      cbc.setLogLevel(loglevel);
      cbc.solver()->messageHandler()->setLogLevel(0);
      clp_->resolve();
      cbc.setStrategy(*strategy_);
      cbc.setLogLevel(loglevel);
      cbc.solver()->messageHandler()->setLogLevel(0);
      cbc.setMaximumSeconds(maxTime);
      cbc.setCutoff(cutoff);
      cbc.setDblParam( CbcModel::CbcAllowableFractionGap, gap_tol_);

      //cbc.solver()->writeMpsNative("FP.mps", NULL, NULL, 1);
      cbc.branchAndBound();
      lowBound_ = cbc.getBestPossibleObjValue();

      if (cbc.isProvenOptimal() || cbc.isProvenInfeasible())
        optimal_ = true;
      else optimal_ = false;

      if (cbc.getSolutionCount()) {
        if (!integerSolution_)
          integerSolution_ = new double[clp_->getNumCols()];
        CoinCopyN(cbc.bestSolution(), clp_->getNumCols(), integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
      nodeCount_ = cbc.getNodeCount();
      iterationCount_ = cbc.getIterationCount();
      delete strat_default;
    }
    else 
#ifdef COIN_HAS_CPX
    if (cpx_) {
      cpx_->switchToMIP();
      CPXENVptr env = cpx_->getEnvironmentPtr();
      CPXLPptr orig_lp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

      int s;
      CPXLPptr cpxlp = CPXcloneprob(env, orig_lp, &s);
      double gap_tol = std::max(0.,gap_tol_- gap_tol_*(1e-01));

#ifdef SHIFT_CUTOFF
      if(cutoff < 1e20){
        cutoff = cutoff-fabs(cutoff)*gap_tol_*0.2;
        gap_tol = gap_tol_*0.8;
      }
#endif

      CPXsetdblparam(env, CPX_PARAM_TILIM, maxTime);
      CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 1);
      CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);
      CPXsetdblparam(env, CPX_PARAM_EPGAP, gap_tol);
      CPXsetintparam( env, CPX_PARAM_PREIND, CPX_ON );

      //CPXwriteprob(env, cpxlp, "OA_trunk","MPS");
      //CPXwriteparam(env, "params.txt");
      

      cpx_->messageHandler()->setLogLevel(loglevel);

      int status = CPXmipopt(env,cpxlp);
      CHECK_CPX_STAT("mipopt",status)

      int stat = CPXgetstat( env, cpxlp);
      bool infeasible = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_ABORT_INFEAS) || (stat == CPXMIP_TIME_LIM_INFEAS) || (stat == CPXMIP_NODE_LIM_INFEAS) || (stat == CPXMIP_FAIL_INFEAS)
                        || (stat == CPXMIP_MEM_LIM_INFEAS) || (stat == CPXMIP_INForUNBD);
      optimal_ = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL) || (stat == CPXMIP_INForUNBD); 
      nodeCount_ = CPXgetnodecnt(env , cpxlp);
      iterationCount_ = CPXgetmipitcnt(env , cpxlp);
      status = CPXgetbestobjval(env, cpxlp, &lowBound_);
      CHECK_CPX_STAT("getbestobjval",status)
       
      if(!infeasible){
         if(!integerSolution_){
           integerSolution_ = new double[cpx_->getNumCols()];
         }
         CPXgetmipx(env, cpxlp, integerSolution_, 0, cpx_->getNumCols() -1);
         CHECK_CPX_STAT("getmipx",status)
      }
      else {
        if (integerSolution_) {
          delete [] integerSolution_;
          integerSolution_ = NULL;
        }
      }
      CPXfreeprob(env, &cpxlp);
      cpx_->switchToLP();
    }
    else {
#else
     {
#endif
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
      }
}

  void
  SubMipSolver::optimize_with_lazy_constraints(double cutoff, int loglevel, double maxTime, const OsiCuts &cs)
  {
    if (clp_) {
      fprintf(stderr, "Function optimize_with_lazy_constraints can only be used with CPLEX\n");
      optimize(cutoff,loglevel, maxTime);
    }
    else 
#ifdef COIN_HAS_CPX
    if (cpx_) {
      cpx_->switchToMIP();
      CPXENVptr env = cpx_->getEnvironmentPtr();
      CPXLPptr cpxlp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

// Remove all the cuts and declare them as lazy constraints

      int orig_nrows = CPXgetnumrows(env, cpxlp) - cs.sizeRowCuts();
      /* printf("Number of rows %i\n", cs.sizeRowCuts()); */
      CPXdelrows(env, cpxlp, orig_nrows, CPXgetnumrows(env, cpxlp) - 1);
      
      int rcnt = cs.sizeRowCuts(), nzcnt = 0;
      vector<double> rhs(rcnt);
      vector<char> sense(rcnt);
      vector<int> beg(rcnt);
      vector<int> ind;
      vector<double> val; 
      double infty = cpx_->getInfinity();
 
      for(int i =0 ; i < rcnt ; i++){
        const OsiRowCut &r = cs.rowCut(i);
        const double lb = r.lb(), ub=r.ub();
        if(ub >= infty) {
          sense[i] = 'G';
          rhs[i] = lb;
        }
        else if (lb <= infty) {
           sense[i] = 'L';
           rhs[i] = ub;
        }
        else {
          assert(lb == ub);
          sense[i] = 'E';
          rhs[i] = ub;
        }
        beg[i] = nzcnt;
        nzcnt += r.row().getNumElements();
      }

      ind.resize(nzcnt);
      val.resize(nzcnt);      
      for(int i =0 ; i < rcnt ; i++){
        const OsiRowCut &r = cs.rowCut(i);
        const double * el = r.row().getElements();
        const int * id = r.row().getIndices();
        int nz = r.row().getNumElements();
        std::copy(el, el + nz, val() + beg[i]);
        std::copy(id, id + nz, ind() + beg[i]);
      }

      CPXaddlazyconstraints(env, cpxlp, rcnt, nzcnt, rhs(), sense(), beg(), ind(), val(), NULL);
      CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);

      CPXsetdblparam(env, CPX_PARAM_TILIM, maxTime);
      CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 1);
      CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);
      CPXsetdblparam(env, CPX_PARAM_EPGAP, gap_tol_);
      cpx_->messageHandler()->setLogLevel(loglevel);
      int status = CPXmipopt(env,cpxlp);
      CHECK_CPX_STAT("mipopt",status)

      int stat = CPXgetstat( env, cpxlp);
      bool infeasible = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_ABORT_INFEAS) || (stat == CPXMIP_TIME_LIM_INFEAS) || (stat == CPXMIP_NODE_LIM_INFEAS) || (stat == CPXMIP_FAIL_INFEAS)
                        || (stat == CPXMIP_MEM_LIM_INFEAS) || (stat == CPXMIP_INForUNBD);
      optimal_ = (stat == CPXMIP_INFEASIBLE) || (stat == CPXMIP_OPTIMAL) || (stat == CPXMIP_OPTIMAL_TOL) || (stat == CPXMIP_INForUNBD);
      nodeCount_ = CPXgetnodecnt(env , cpxlp);
      iterationCount_ = CPXgetmipitcnt(env , cpxlp);
      status = CPXgetbestobjval(env, cpxlp, &lowBound_);
      CHECK_CPX_STAT("getbestobjval",status)
       
      if(!infeasible){
         if(!integerSolution_){
           integerSolution_ = new double[cpx_->getNumCols()];
         }
         CPXgetmipx(env, cpxlp, integerSolution_, 0, cpx_->getNumCols() -1);
         CHECK_CPX_STAT("getmipx",status)
      }
      else {
        if (integerSolution_) {
          delete [] integerSolution_;
          integerSolution_ = NULL;
        }
      }
      cpx_->switchToLP();
      CPXfreelazyconstraints(env, cpxlp);
      CPXaddrows(env, cpxlp, 0, rcnt, nzcnt, rhs(), sense(), beg(), ind(), val(), NULL, NULL);
    }
    else {
#else
     {
#endif
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
      }
}

   /** Assign a strategy. */
   void 
   SubMipSolver::setStrategy(CbcStrategyDefault * strategy)
   {
     if (strategy_) delete strategy_;
     strategy_ = dynamic_cast<CbcStrategyDefault *>(strategy->clone());
     assert(strategy_);
   }

  /** Register options.*/
  void
  SubMipSolver::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("MILP Solver", RegisteredOptions::BonminCategory);
    roptions->AddStringOption3("milp_solver",
        "Choose the subsolver to solve MILP sub-problems in OA decompositions.",
        "Cbc_D",
        "Cbc_D","Coin Branch and Cut with its default",
        "Cbc_Par", "Coin Branch and Cut with passed parameters",
        "Cplex","IBM Cplex",
        " To use Cplex, a valid license is required and you should have compiled OsiCpx in COIN-OR  (see Osi documentation).");
    roptions->setOptionExtraInfo("milp_solver",64);

    roptions->AddBoundedIntegerOption("cpx_parallel_strategy",
                           "Strategy of parallel search mode in CPLEX.",
                           -1, 1, 0,
                           "-1 = opportunistic, 0 = automatic, 1 = deterministic (refer to CPLEX documentation)"
                           );
    roptions->setOptionExtraInfo("cpx_parallel_strategy",64);

    roptions->AddLowerBoundedIntegerOption("number_cpx_threads",
                           "Set number of threads to use with cplex.",
                           0, 0,
                           "(refer to CPLEX documentation)"
                           );
    roptions->setOptionExtraInfo("number_cpx_threads",64);

    
    roptions->AddStringOption2("milp_strategy",
        "Choose a strategy for MILPs.",
        "solve_to_optimality",
        "find_good_sol","Stop sub milps when a solution improving the incumbent is found",
        "solve_to_optimality", "Solve MILPs to optimality",
        "");
    roptions->setOptionExtraInfo("milp_strategy",64);

    roptions->SetRegisteringCategory("Output and Loglevel", RegisteredOptions::BonminCategory);
    roptions->AddBoundedIntegerOption("milp_log_level",
        "specify MILP solver log level.",
        0,4,0,
        "Set the level of output of the MILP subsolver in OA : "
        "0 - none, 1 - minimal, 2 - normal low, 3 - normal high"
                                     );
    roptions->setOptionExtraInfo("milp_log_level",64);

  }
}/* Ends Bonmin namespace.*/
