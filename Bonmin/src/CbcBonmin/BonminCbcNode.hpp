// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// John J. Forrest, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#ifndef BonminCbcNode_H
#define BonminCbcNode_H

#include "CbcNode.hpp"



/** \brief Holds information for recreating a subproblem by incremental change
	   from the parent for Bonmin

  A BonminBonminCbcPartialNodeInfo object contains changes to the bounds and basis, and
  additional cuts, required to recreate a subproblem by modifying and
  augmenting the parent subproblem.
*/

class BonminCbcFullNodeInfo : public CbcFullNodeInfo
{

public:
  friend class BonminCbcPartialNodeInfo;
  // Default Constructor
  BonminCbcFullNodeInfo ();

  // Constructor from current state
  BonminCbcFullNodeInfo (CbcModel * model, int numberRowsAtContinuous);

  // Copy constructor
  BonminCbcFullNodeInfo ( const BonminCbcFullNodeInfo &);

  // Destructor
  ~BonminCbcFullNodeInfo ();

  /// Clone
  virtual CbcNodeInfo * clone() const;

  /**Method called when all direct sons have been explored to flush
     useless warm start information.*/
  virtual void allBranchesGone();

  /** Number of consecutive infeasible parents only recorded if node is infeasible*/
  inline int getSequenceOfInfeasiblesSize()
  {
    return sequenceOfInfeasiblesSize_;
  }
  /** Number of consecutive unsolved parents only recorded if node is infeasible*/
  inline int getSequenceOfUnsolvedSize()
  {
    return sequenceOfUnsolvedSize_;
  }
private:
  /* Data values */
  /** Number of consecutive infeasible parents only recorded if node is infeasible*/
  int sequenceOfInfeasiblesSize_;
  /** Number of consecutive unsolved parents only recorded if node is infeasible*/
  int sequenceOfUnsolvedSize_;
private:

  /// Illegal Assignment operator
  BonminCbcFullNodeInfo & operator=(const BonminCbcFullNodeInfo& rhs);
};

/** \brief Holds information for recreating a subproblem by incremental change
	   from the parent for Bonmin

  A BonminBonminCbcPartialNodeInfo object contains changes to the bounds and basis, and
  additional cuts, required to recreate a subproblem by modifying and
  augmenting the parent subproblem.
*/

class BonminCbcPartialNodeInfo : public CbcPartialNodeInfo
{

public:
  // Default Constructor
  BonminCbcPartialNodeInfo ();

  // Constructor from current state
  BonminCbcPartialNodeInfo (CbcModel * model, CbcNodeInfo * parent, CbcNode * owner,
      int numberChangedBounds,const int * variables,
      const double * boundChanges,
      const CoinWarmStartDiff *basisDiff) ;

  // Copy constructor
  BonminCbcPartialNodeInfo ( const BonminCbcPartialNodeInfo &);

  // Destructor
  ~BonminCbcPartialNodeInfo ();

  /// Clone
  virtual CbcNodeInfo * clone() const;

  /**Method called when all direct sons have been explored to flush
     useless warm start information.*/
  virtual void allBranchesGone();

  /** Number of consecutive infeasible parents only recorded if node is infeasible*/
  inline int getSequenceOfInfeasiblesSize()
  {
    return sequenceOfInfeasiblesSize_;
  }
  /** Number of consecutive unsolved parents only recorded if node is infeasible*/
  inline int getSequenceOfUnsolvedSize()
  {
    return sequenceOfUnsolvedSize_;
  }
private:
  /* Data values */
  /** Number of consecutive infeasible parents only recorded if node is infeasible*/
  int sequenceOfInfeasiblesSize_;
  /** Number of consecutive unsolved parents only recorded if node is infeasible*/
  int sequenceOfUnsolvedSize_;
private:

  /// Illegal Assignment operator
  BonminCbcPartialNodeInfo & operator=(const BonminCbcPartialNodeInfo& rhs);
};
#endif
