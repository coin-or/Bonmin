// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// John J. Forrest, International Business Machines Corporation
// P. Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "BonminCbcNlpStrategy.hpp"
#include "BonminCbcNode.hpp"

#include "IpoptInterface.hpp"
namespace Bonmin{
// Default Constructor
BonminCbcNlpStrategy::BonminCbcNlpStrategy(int maxFailures,
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
BonminCbcNlpStrategy::~BonminCbcNlpStrategy ()
{}

// Clone
CbcStrategy *
BonminCbcNlpStrategy::clone() const
{
  return new BonminCbcNlpStrategy(*this);
}

// Copy constructor
BonminCbcNlpStrategy::BonminCbcNlpStrategy(const BonminCbcNlpStrategy & rhs)
    :
    hasFailed_(false),
    maxFailure_(rhs.maxFailure_),
    maxInfeasible_(rhs.maxInfeasible_),
    pretendFailIsInfeasible_(rhs.pretendFailIsInfeasible_)
{}
// Return a new Full node information pointer (descendant of CbcFullNodeInfo)
CbcNodeInfo *
BonminCbcNlpStrategy::fullNodeInfo(CbcModel * model,int numberRowsAtContinuous) const
{
  return new CbcFullNodeInfo(model,numberRowsAtContinuous);
}
// Return a new Partial node information pointer (descendant of CbcPartialNodeInfo)
CbcNodeInfo *
BonminCbcNlpStrategy::partialNodeInfo(CbcModel * model, CbcNodeInfo * parent, CbcNode * owner,
    int numberChangedBounds,const int * variables,
    const double * boundChanges,
    const CoinWarmStartDiff *basisDiff) const
{
  return new BonminCbcPartialNodeInfo(model,parent, owner, numberChangedBounds, variables,
      boundChanges,basisDiff);
}
/* After a CbcModel::resolve this can return a status
   -1 no effect
   0 treat as optimal
   1 as 0 but do not do any more resolves (i.e. no more cuts)
   2 treat as infeasible
*/
int
BonminCbcNlpStrategy::status(CbcModel * model, CbcNodeInfo * parent,int whereFrom)
{
  OsiSolverInterface * solver = model->solver();//get solver
  int feasible = 1;
  bool solved = true;
  int returnStatus = -1;
  BonminCbcPartialNodeInfo * bmNodeInfo = dynamic_cast<BonminCbcPartialNodeInfo *>(parent);
  if(!bmNodeInfo) return -1;

  int seqOfInfeasiblesSize = bmNodeInfo->getSequenceOfInfeasiblesSize();
  int seqOfUnsolvedSize = bmNodeInfo->getSequenceOfUnsolvedSize();


  if(solver->isAbandoned()) {
    solved = false;
    seqOfUnsolvedSize++;
    ;
  }
  else if(solver->isProvenPrimalInfeasible()) {
    feasible = 0;
    seqOfInfeasiblesSize++;
  }

  if((seqOfUnsolvedSize==0) || (maxFailure_ == 0) &&
      (maxInfeasible_== 0) || (seqOfInfeasiblesSize==0))

    if(feasible && seqOfInfeasiblesSize > 1) {
      std::cerr<<"Feasible node while father was infeasible."
      <<std::endl;
    }

  if(solved && seqOfUnsolvedSize > 1) {
    std::cerr<<"Solved node while father was unsolved."
    <<std::endl;
  }

  if(seqOfInfeasiblesSize < maxInfeasible_ &&
      solved && !feasible) {
    std::cerr<<"Branching on infeasible node, sequence of infeasibles size "
    <<seqOfInfeasiblesSize<<std::endl;
    // Have to make sure that we will branch
    IpoptInterface * ipopt = dynamic_cast<IpoptInterface *>(solver);
    ipopt->forceBranchable();
    //change objective value
    returnStatus = 0;

  }

  if(!solved && parent != NULL &&
      seqOfUnsolvedSize <= maxFailure_) {
    std::cout<<"Branching on unsolved node, sequence of unsolved size "<<seqOfUnsolvedSize<<std::endl;
    // Have to make sure that we will branch
    IpoptInterface * ipopt = dynamic_cast<IpoptInterface *>(solver);
    ipopt->forceBranchable();     //      feasible=1;
    returnStatus = 0;
  }

  if(solver->isAbandoned() && parent != NULL &&
      seqOfUnsolvedSize > maxFailure_) {
    hasFailed_ = true;
    IpoptInterface * ipopt = dynamic_cast<IpoptInterface *>(solver);
    if(pretendFailIsInfeasible_) {
      //force infeasible
      ipopt->forceInfeasible();
      returnStatus = 2;
    }
    else
      throw ipopt->newUnsolvedError(ipopt->getOptStatus());
  }
  return returnStatus;
}

void
BonminCbcNlpStrategy::setupCutGenerators(CbcModel &model)
{}

void
BonminCbcNlpStrategy::setupHeuristics(CbcModel &model)
{}

void
BonminCbcNlpStrategy::setupPrinting(CbcModel &model, int toto)
{}

void
BonminCbcNlpStrategy::setupOther(CbcModel &model)
{}
}
