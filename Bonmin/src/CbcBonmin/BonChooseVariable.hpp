// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonChooseVariable_H
#define BonChooseVariable_H

#include "OsiChooseVariable.hpp"
#include "BonCurvatureEstimator.hpp"
#include "BonOsiTMINLPInterface.hpp"

namespace Bonmin {

/** This class chooses a variable to branch on

    The base class just chooses the variable and direction without strong branching but it 
    has information which would normally be used by strong branching e.g. to re-enter
    having fixed a variable but using same candidates for strong branching.

    The flow is :
    a) initialize the process.  This decides on strong branching list
       and stores indices of all infeasible objects  
    b) do strong branching on list.  If list is empty then just
       choose one candidate and return without strong branching.  If not empty then
       go through list and return best.  However we may find that the node is infeasible
       or that we can fix a variable.  If so we return and it is up to user to call
       again (after fixing a variable).
*/

class BonChooseVariable : public OsiChooseVariable  {
 
public:

  /// Constructor from solver (so we can set up arrays etc)
  BonChooseVariable (OsiTMINLPInterface * solver);

  /// Copy constructor 
  BonChooseVariable (const BonChooseVariable &);
   
  /// Assignment operator 
  BonChooseVariable & operator= (const BonChooseVariable& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~BonChooseVariable ();

  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. */
  virtual int setupList ( OsiBranchingInformation *info, bool initialize);
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

#ifdef AWAWAWAW
  /// Given a candidate fill in useful information e.g. estimates
  virtual void updateInformation( const OsiBranchingInformation *info,
				  int branch, OsiHotInfo * hotInfo);
#endif

protected:
  SmartPtr<CurvatureEstimator> cur_estimator_;

private:
  /// Default Constructor 
  BonChooseVariable ();

};

}

#endif
