// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonQPStrongBranching_H
#define BonQPStrongBranching_H

#include "OsiChooseVariable.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonBranchingTQP.hpp"

namespace Bonmin {

/** This class chooses a variable to branch on

    This implementation solves the QP model for different branches
    (strong branching).
*/

class BonQPStrongBranching : public OsiChooseVariable  {
 
public:

  /// Constructor from solver (so we can set up arrays etc)
  BonQPStrongBranching (OsiTMINLPInterface * solver);

  /// Copy constructor 
  BonQPStrongBranching (const BonQPStrongBranching &);
   
  /// Assignment operator 
  BonQPStrongBranching & operator= (const BonQPStrongBranching& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~BonQPStrongBranching ();

#define UseOurOwn
#ifdef UseOurOwn
  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. */
  virtual int setupList ( OsiBranchingInformation *info, bool initialize);
#endif
  /** Choose a variable
      Returns - 
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from bestObjectIndex() and bestWhichWay()
     We can pick up a forced branch (can change bound) from firstForcedObjectIndex() and firstForcedWhichWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
     If fixVariables is true then 2,3,4 are all really same as problem changed
  */
  virtual int chooseVariable( OsiSolverInterface * solver, OsiBranchingInformation *info, bool fixVariables);
  /**  This is a utility function which does strong branching on
       a list of objects and stores the results in OsiHotInfo.objects.
       On entry the object sequence is stored in the OsiHotInfo object
       and maybe more.
       It returns -
       -1 - one branch was infeasible both ways
       0 - all inspected - nothing can be fixed
       1 - all inspected - some can be fixed (returnCriterion==0)
       2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
       3 - returning because max time
       
  */

protected:
  SmartPtr<BranchingTQP> branching_tqp_;

private:
  /// Default Constructor 
  BonQPStrongBranching ();

};

}

#endif
