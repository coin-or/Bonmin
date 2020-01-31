// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// John J. Forrest, International Business Machines Corporation
// P. Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006


#include <string>
#include <cassert>
#include <cfloat>
#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CbcModel.hpp"
#include "BonCbcNode.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptWarmStart.hpp"
#include "BonOsiTMINLPInterface.hpp"

using namespace std;


namespace Bonmin
{
//Default constructor
  BonCbcFullNodeInfo::BonCbcFullNodeInfo()
      :
      CbcFullNodeInfo(),
      sequenceOfInfeasiblesSize_(0),
      sequenceOfUnsolvedSize_(0)
  {}

  BonCbcFullNodeInfo::BonCbcFullNodeInfo(CbcModel * model,
      int numberRowsAtContinuous) :
      CbcFullNodeInfo(model, numberRowsAtContinuous),
      sequenceOfInfeasiblesSize_(0),
      sequenceOfUnsolvedSize_(0)
  {}

// Copy constructor
  BonCbcFullNodeInfo::BonCbcFullNodeInfo ( const BonCbcFullNodeInfo &other):
      CbcFullNodeInfo(other),
      sequenceOfInfeasiblesSize_(other.sequenceOfInfeasiblesSize_),
      sequenceOfUnsolvedSize_(other.sequenceOfUnsolvedSize_)

  {}


  void
  BonCbcFullNodeInfo::allBranchesGone()
  {
    IpoptWarmStart * ipws = dynamic_cast<IpoptWarmStart *>(basis_);
    if (ipws)
      ipws->flushPoint();
  }

  BonCbcFullNodeInfo::~BonCbcFullNodeInfo()
{}

  CbcNodeInfo *
  BonCbcFullNodeInfo::clone() const
  {
    return new BonCbcFullNodeInfo(*this);
  }

  void
  BonCbcFullNodeInfo::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {


    roptions->SetRegisteringCategory("Nonconvex problems", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedIntegerOption("max_consecutive_infeasible",
        "Number of consecutive infeasible subproblems before aborting a"
        " branch.",
        0,0,
        "Will continue exploring a branch of the tree until \"max_consecutive_infeasible\""
        "consecutive problems are locally infeasible by the NLP sub-solver.");
    roptions->setOptionExtraInfo("max_consecutive_infeasible",8);

    roptions->SetRegisteringCategory("NLP solution robustness", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedIntegerOption
    ("max_consecutive_failures",
     "(temporarily removed) Number $n$ of consecutive unsolved problems before aborting a branch of the tree.",
     0,10,
     "When $n > 0$, continue exploring a branch of the tree until $n$ "
     "consecutive problems in the branch are unsolved (we call unsolved a problem for which Ipopt can not "
     "guarantee optimality within the specified tolerances).");
    roptions->setOptionExtraInfo("max_consecutive_failures",8);

  }


  /****************************************************************************************************/

// Default constructor
  BonCbcPartialNodeInfo::BonCbcPartialNodeInfo ()
      : CbcPartialNodeInfo(),
      sequenceOfInfeasiblesSize_(0),
      sequenceOfUnsolvedSize_(0)
  {}
// Constructor from current state
  BonCbcPartialNodeInfo::BonCbcPartialNodeInfo (CbcModel * model,CbcNodeInfo *parent, CbcNode *owner,
      int numberChangedBounds,
      const int *variables,
      const double *boundChanges,
      const CoinWarmStartDiff *basisDiff)
      : CbcPartialNodeInfo(parent,owner,numberChangedBounds,variables,
          boundChanges,basisDiff),
      sequenceOfInfeasiblesSize_(0),
      sequenceOfUnsolvedSize_(0)
  {
    BonCbcPartialNodeInfo * nlpParent = dynamic_cast<BonCbcPartialNodeInfo *> (parent);
    int numberInfeasible = 0;
    int numberUnsolved = 0;
    if (nlpParent)//father is not root
    {
      numberInfeasible = nlpParent->getSequenceOfInfeasiblesSize();
      numberUnsolved =  nlpParent->getSequenceOfUnsolvedSize();
//       if(!nlpParent->numberBranchesLeft_){
// 	IpoptWarmStartDiff * ipws = dynamic_cast<IpoptWarmStartDiff *>(nlpParent->basisDiff_);
// 	ipws->flushPoint();
//       }
    }
    else {
      BonCbcFullNodeInfo * nlpRoot = dynamic_cast<BonCbcFullNodeInfo *> (parent);
      if (nlpRoot) {
        numberInfeasible = nlpRoot->getSequenceOfInfeasiblesSize();
        numberUnsolved =  nlpRoot->getSequenceOfUnsolvedSize();
      }
    }
    if (model->solver()->isAbandoned() ||
        model->solver()->isIterationLimitReached())
      sequenceOfUnsolvedSize_ = numberUnsolved + 1;

    if (model->solver()->isProvenPrimalInfeasible())
      sequenceOfInfeasiblesSize_ = numberInfeasible + 1;
  }

  BonCbcPartialNodeInfo::BonCbcPartialNodeInfo (const BonCbcPartialNodeInfo & rhs)

      : CbcPartialNodeInfo(rhs),
      sequenceOfInfeasiblesSize_(rhs.sequenceOfInfeasiblesSize_),
      sequenceOfUnsolvedSize_(rhs.sequenceOfUnsolvedSize_)

{}

  CbcNodeInfo *
  BonCbcPartialNodeInfo::clone() const
  {
    return (new BonCbcPartialNodeInfo(*this));
  }

  void
  BonCbcPartialNodeInfo::allBranchesGone()
  {
    IpoptWarmStartDiff * ipws = dynamic_cast<IpoptWarmStartDiff *>(basisDiff_);
    if (ipws)
      ipws->flushPoint();
  }

  BonCbcPartialNodeInfo::~BonCbcPartialNodeInfo ()
{}
}
