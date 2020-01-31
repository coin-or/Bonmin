// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// John J. Forrest, International Business Machines Corporation
// P. Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#include <cassert>
#include <cmath>
#include <cfloat>

#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "BonCbcNlpStrategy.hpp"
#include "BonCbcNode.hpp"

#include "BonOsiTMINLPInterface.hpp"
namespace Bonmin
{
// Default Constructor
  CbcNlpStrategy::CbcNlpStrategy(int maxFailures,
      int maxInfeasibles,
      int pretendFailIsInfeasible)
      :
      hasFailed_(false),
      maxFailure_(maxFailures),
      maxInfeasible_(maxInfeasibles),
      pretendFailIsInfeasible_(pretendFailIsInfeasible)
  {
    setPreProcessState(0);
  }


// Destructor
  CbcNlpStrategy::~CbcNlpStrategy ()
  {}

// Clone
  CbcStrategy *
  CbcNlpStrategy::clone() const
  {
    return new CbcNlpStrategy(*this);
  }

// Copy constructor
  CbcNlpStrategy::CbcNlpStrategy(const CbcNlpStrategy & rhs)
      :
      hasFailed_(false),
      maxFailure_(rhs.maxFailure_),
      maxInfeasible_(rhs.maxInfeasible_),
      pretendFailIsInfeasible_(rhs.pretendFailIsInfeasible_)
  {}
// Return a new Full node information pointer (descendant of CbcFullNodeInfo)
  CbcNodeInfo *
  CbcNlpStrategy::fullNodeInfo(CbcModel * model,int numberRowsAtContinuous) const
  {
    return new CbcFullNodeInfo(model,numberRowsAtContinuous);
  }
// Return a new Partial node information pointer (descendant of CbcPartialNodeInfo)
  CbcNodeInfo *
  CbcNlpStrategy::partialNodeInfo(CbcModel * model, CbcNodeInfo * parent, CbcNode * owner,
      int numberChangedBounds,const int * variables,
      const double * boundChanges,
      const CoinWarmStartDiff *basisDiff) const
  {
    return new BonCbcPartialNodeInfo(model,parent, owner, numberChangedBounds, variables,
        boundChanges,basisDiff);
  }
  /* After a CbcModel::resolve this can return a status
     -1 no effect
     0 treat as optimal
     1 as 0 but do not do any more resolves (i.e. no more cuts)
     2 treat as infeasible
  */
  int
  CbcNlpStrategy::status(CbcModel * model, CbcNodeInfo * parent,int whereFrom)
  {

    OsiSolverInterface * solver = model->solver();//get solver
    int feasible = 1;
    bool solved = true;
    int returnStatus = -1;
    BonCbcPartialNodeInfo * bmNodeInfo = dynamic_cast<BonCbcPartialNodeInfo *>(parent);
    if (!bmNodeInfo) return -1;

    int seqOfInfeasiblesSize = bmNodeInfo->getSequenceOfInfeasiblesSize();
    int seqOfUnsolvedSize = bmNodeInfo->getSequenceOfUnsolvedSize();


    if (solver->isAbandoned()) {
      solved = false;
      seqOfUnsolvedSize++;
      ;
    }
    else if (solver->isProvenPrimalInfeasible()) {
      feasible = 0;
      seqOfInfeasiblesSize++;
    }

    if (((seqOfUnsolvedSize==0) || (maxFailure_ == 0)) &&
        ((maxInfeasible_== 0) || (seqOfInfeasiblesSize==0)))

      if (feasible && seqOfInfeasiblesSize > 1) {
        (*model->messageHandler())<<"Feasible node while father was infeasible."
        <<CoinMessageEol;
      }

    if (solved && seqOfUnsolvedSize > 1) {
      (*model->messageHandler())<<"Solved node while father was unsolved."
      <<CoinMessageEol;
    }

    if (seqOfInfeasiblesSize < maxInfeasible_ &&
        solved && !feasible) {
      (*model->messageHandler())<<"Branching on infeasible node, sequence of infeasible size "
      <<seqOfInfeasiblesSize<<CoinMessageEol;
      // Have to make sure that we will branch
      OsiTMINLPInterface * ipopt = dynamic_cast<OsiTMINLPInterface *>(solver);
      ipopt->forceBranchable();
      //change objective value
      returnStatus = 0;

    }

    if (!solved && parent != NULL &&
        seqOfUnsolvedSize <= maxFailure_) {
      (*model->messageHandler())<<"Branching on unsolved node, sequence of unsolved size "<<seqOfUnsolvedSize<<CoinMessageEol;
      // Have to make sure that we will branch
      OsiTMINLPInterface * osiMinlp = dynamic_cast<OsiTMINLPInterface *>(solver);
      osiMinlp->forceBranchable();     //      feasible=1;
      returnStatus = 0;
    }

    if (solver->isAbandoned() && parent != NULL &&
        seqOfUnsolvedSize > maxFailure_) {
      hasFailed_ = true;
      OsiTMINLPInterface * osiMinlp =
        dynamic_cast<OsiTMINLPInterface *>(solver);
      if (pretendFailIsInfeasible_) {
        //force infeasible
        osiMinlp->forceInfeasible();
        returnStatus = 2;
      }
      else {
        std::string probName;
        osiMinlp->getStrParam(OsiProbName,probName);
        throw osiMinlp->newUnsolvedError(0, osiMinlp->problem(), probName);
      }
    }
    return returnStatus;
  }

  void
  CbcNlpStrategy::setupCutGenerators(CbcModel &model)
  {}

  void
  CbcNlpStrategy::setupHeuristics(CbcModel &model)
  {}

  void
  CbcNlpStrategy::setupPrinting(CbcModel &model, int toto)
  {}

  void
  CbcNlpStrategy::setupOther(CbcModel &model)
  {}
}
