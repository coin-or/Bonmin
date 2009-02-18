// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006


// Code separated from BonOaDecBase to try to clarify OAs
#include "BonSubMipSolver.hpp"
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiClpSolverInterface.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif

namespace Bonmin {
  /** Constructor */
  SubMipSolver::SubMipSolver(OsiSolverInterface * lp,
      const CbcStrategy * strategy):
      lp_(lp),
      clp_(NULL),
      cpx_(NULL),
      cbc_(NULL),
      lowBound_(-COIN_DBL_MAX),
      optimal_(false),
      integerSolution_(NULL),
      strategy_(NULL)
  {
    clp_ = (lp_ == NULL)? NULL :
        dynamic_cast<OsiClpSolverInterface *>(lp_);
#ifdef COIN_HAS_CPX
    cpx_ = (lp_ == NULL)? NULL :
        dynamic_cast<OsiCpxSolverInterface *>(lp_);
#endif
    if (strategy) strategy_ = strategy->clone();
  }
  SubMipSolver::~SubMipSolver()
  {
    if (strategy_) delete strategy_;
    if (integerSolution_) delete [] integerSolution_;
    if (cbc_) delete cbc_;
  }

  /** Assign lp solver. */
  void
  SubMipSolver::setLpSolver(OsiSolverInterface * lp)
  {
    lp_ = lp;
    clp_ = (lp_ == NULL) ? NULL :
        dynamic_cast<OsiClpSolverInterface *>(lp_);
#ifdef COIN_HAS_CPX
    cpx_ = (lp_ == NULL) ? NULL :
        dynamic_cast<OsiCpxSolverInterface *>(lp_);
#endif
    lowBound_ = -COIN_DBL_MAX;
    optimal_ = false;
    if (integerSolution_) {
      delete [] integerSolution_;
      integerSolution_ = NULL;
    }
  }


 void 
 SubMipSolver::find_good_sol(double cutoff, int loglevel, double max_time,
                                    int max_node){

     if(clp_){
      if (!strategy_)
        strategy_ = new CbcStrategyDefault(1,0,0, loglevel);

      OsiBabSolver empty;
      if (cbc_) delete cbc_;
      cbc_ = new CbcModel(*clp_);
      cbc_->solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(cbc_->messagesPointer()->source_,"OaCbc");

      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      clp_->resolve();
      cbc_->setStrategy(*strategy_);
      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      cbc_->setMaximumNodes(max_node);
      cbc_->setMaximumSeconds(max_time);
      cbc_->setMaximumSolutions(1);
      cbc_->setCutoff(cutoff);

      cbc_->branchAndBound();
      lowBound_ = cbc_->getBestPossibleObjValue();

      if (cbc_->isProvenOptimal() || cbc_->isProvenInfeasible())
        optimal_ = true;
      else optimal_ = false;

      if (cbc_->getSolutionCount()) {
        if (!integerSolution_)
          integerSolution_ = new double[lp_->getNumCols()];
        CoinCopyN(cbc_->bestSolution(), lp_->getNumCols(), integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
      nodeCount_ = cbc_->getNodeCount();
      iterationCount_ = cbc_->getIterationCount();

     }
     else if (cpx_){
#ifndef COIN_HAS_CPX
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
#else
        CPXENVptr env = cpx_->getEnvironmentPtr();
        CPXLPptr cpxlp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
        CPXsetintparam(env, CPX_PARAM_NODELIM, max_node);
        CPXsetdblparam(env, CPX_PARAM_TILIM, max_time);
        CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);

        CPXsetintparam(env,CPX_PARAM_INTSOLLIM, 10000);
        CPXsetintparam(env,CPX_PARAM_NODELIM, 100000);
        CPXsetdblparam(env,CPX_PARAM_TILIM, max_time);

        cpx_->branchAndBound();

        int status = CPXgetbestobjval(env, cpxlp, &lowBound_);
     
        int stat = CPXgetstat( env, cpxlp);
        optimal_ |= (stat == CPXMIP_INFEASIBLE); 
        nodeCount_ = CPXgetnodecnt(env , cpxlp);
        iterationCount_ = CPXgetmipitcnt(env , cpxlp);

        CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
        CPXsetintparam(env,CPX_PARAM_NODELIM, 1000);
        while(stat == CPXMIP_NODE_LIM_INFEAS){
           cpx_->branchAndBound();
           int stat = CPXgetstat( env, cpxlp);
           optimal_ |= (stat == CPXMIP_INFEASIBLE);
           nodeCount_ = CPXgetnodecnt(env , cpxlp);
           iterationCount_ = CPXgetmipitcnt(env , cpxlp);
        }
      
        if (status)
          throw CoinError("Error in getting some CPLEX information",
                          "OaDecompositionBase::SubMipSolver",
                          "performLocalSearch");
#endif
    }
  }

  void
  SubMipSolver::optimize(double cutoff, int loglevel, double maxTime,
      int maxNodes)
  {
    if (clp_) {
      if (!strategy_)
        strategy_ = new CbcStrategyDefault(1,0,0, loglevel);

      OsiBabSolver empty;
      if (cbc_) delete cbc_;
      cbc_ = new CbcModel(*clp_);
      cbc_->solver()->setAuxiliaryInfo(&empty);

      //Change Cbc messages prefixes
      strcpy(cbc_->messagesPointer()->source_,"OaCbc");

      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      clp_->resolve();
      cbc_->setStrategy(*strategy_);
      cbc_->setLogLevel(loglevel);
      cbc_->solver()->messageHandler()->setLogLevel(0);
      cbc_->setMaximumNodes(maxNodes);
      cbc_->setMaximumSeconds(maxTime);
      cbc_->setCutoff(cutoff);

      cbc_->branchAndBound();
      lowBound_ = cbc_->getBestPossibleObjValue();

      if (cbc_->isProvenOptimal() || cbc_->isProvenInfeasible())
        optimal_ = true;
      else optimal_ = false;

      if (cbc_->getSolutionCount()) {
        if (!integerSolution_)
          integerSolution_ = new double[lp_->getNumCols()];
        CoinCopyN(cbc_->bestSolution(), lp_->getNumCols(), integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
      nodeCount_ = cbc_->getNodeCount();
      iterationCount_ = cbc_->getIterationCount();
    }
    else {
      lp_->messageHandler()->setLogLevel(loglevel);
#ifdef COIN_HAS_CPX
      if (cpx_) {
        CPXENVptr env = cpx_->getEnvironmentPtr();
        CPXsetintparam(env, CPX_PARAM_NODELIM, maxNodes);
        CPXsetdblparam(env, CPX_PARAM_TILIM, maxTime);
        CPXsetdblparam(env, CPX_PARAM_CUTUP, cutoff);
        //CpxModel = cpx_;
      }
      else
#endif 
     {
        throw CoinError("Unsuported solver, for local searches you should use clp or cplex",
            "performLocalSearch",
            "OaDecompositionBase::SubMipSolver");
    }

    lp_->branchAndBound();

    optimal_ = lp_->isProvenOptimal();
#ifdef COIN_HAS_CPX
    if (cpx_) {
      //CpxModel = NULL;
      CPXENVptr env = cpx_->getEnvironmentPtr();
      CPXLPptr cpxlp = cpx_->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);

      int status = CPXgetbestobjval(env, cpxlp, &lowBound_);
     
      int stat = CPXgetstat( env, cpxlp);
      optimal_ |= (stat == CPXMIP_INFEASIBLE); 
       nodeCount_ = CPXgetnodecnt(env , cpxlp);
      iterationCount_ = CPXgetmipitcnt(env , cpxlp);
      if (status)
        throw CoinError("Error in getting some CPLEX information","OaDecompositionBase::SubMipSolver","performLocalSearch");
    }
#endif

      if (lp_->getFractionalIndices().size() == 0) {
        if (!integerSolution_)
          integerSolution_ = new double[lp_->getNumCols()];
        CoinCopyN(lp_->getColSolution(), lp_->getNumCols() , integerSolution_);
      }
      else if (integerSolution_) {
        delete [] integerSolution_;
        integerSolution_ = NULL;
      }
    }
  }

   /** Assign a strategy. */
   void 
   SubMipSolver::setStrategy(CbcStrategy * strategy)
   {
     if (strategy_) delete strategy_;
     strategy_ = strategy->clone();
   }

}/* Ends Bonmin namespace.*/
