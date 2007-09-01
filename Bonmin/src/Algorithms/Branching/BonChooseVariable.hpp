// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonChooseVariable_H
#define BonChooseVariable_H

#include "OsiChooseVariable.hpp"
#include "BonCurvatureEstimator.hpp"
#include "BonOsiTMINLPInterface.hpp"

// Forward declaration
class CbcModel;

namespace Bonmin {

  class GuessHeuristic;

/** This class chooses a variable to branch on

    This is the base class for the branching rules in Bonmin (inherits
    from OsiChooseVariable). This class implements a simple strong branching algorithm where the changes in the objective
    value induced by branching on a specific object are estimated with the pure virtual function fill_changes.
*/

class BonChooseVariable : public OsiChooseStrong  {
  /** Statuses for strong branching candidates.*/
  enum StrongStatus{
    NotDone=-1,
    Feasible/** Child is proven feasible.*/,
    Infeasible /** Child is proven infeasible.*/,
    NotFinished /** Child is not finished.*/};
public:
  enum DoStrongReturnStatuses{
  provenInfeasible = -1 /** One branch has two infeasible childs.*/,
  doneNoFixing /** All done no variable can be fixed.*/,
  doneCanFix /** Several variable can be fixed.*/,
  interuptedCanFix /** Interupted and found a variable to fix.*/,
  maxTime /** Interupted because of time limit.*/};

  /** Return statuses for chooseVariable.*/
  enum chooseVariableReturnStatuses{
    infeasibleNode = -1/** Node has been proven infeasible.*/,
    hasCandidate /** Normal termination, found a variable to branch on.*/,
    feasibleNode /** All variable are feasible, the node is feasible.*/,
    canFixAndStrongBranch /** Found variable to fix and also has remaining candidate for strong branching.*/,
    canFixAndBranch/** Found variable to fix and also has a (non-strong) branching candidate.*/,
    canFixNoCandidate /** Can fix variables but does not have strong branching candidates.*/
  };
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

  /// Given a candidate fill in useful information e.g. estimates
  virtual void updateInformation(const OsiBranchingInformation *info,
				  int branch, OsiHotInfo * hotInfo);
  /// Given a branch fill in useful information e.g. estimates
  virtual void updateInformation( int whichObject, int branch, 
				  double changeInObjective, double changeInValue,
				  int status);

  /** Method for setting CbcModel, which is used to get statusOfSearch */
  void setCbcModel(CbcModel* cbc_model)
  {
    cbc_model_ = cbc_model;
  }

  void setOnlyPseudoWhenTrusted(bool only_pseudo_when_trusted)
  {
    only_pseudo_when_trusted_ = only_pseudo_when_trusted;
  }

  /** @name Accessor methods to pseudo cost data*/
  //@{
  int numberObjects() const {
    return numberObjects_;
  }
  const double* upTotalChange() const
  {
    return upTotalChange_;
  }
  const double* downTotalChange() const
  {
    return downTotalChange_;
  }
  const int* upNumber() const 
  {
    return upNumber_;
  }
  const int* downNumber() const 
  {
    return downNumber_;
  }
  //@}

  /// For now, we need to communicate pseudo costs to
  /// GuessHeuristic. Right now, we need to call it back so that it
  /// gets the correct pointers to the arrays (after cloning)
  void registerGuessHeuristic(GuessHeuristic* guessHeuristic)
  {
    guessHeuristic_ = guessHeuristic;
  }

protected:

  /// Holding on the a pointer to the journalist
  SmartPtr<Journalist> jnlst_;

  /// verbosity level
  int bb_log_level_;

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
  virtual int doBonStrongBranching( OsiSolverInterface * solver, 
				    OsiBranchingInformation *info, int numberToDo,
				    OsiHotInfo * results, int returnCriterion,
				    int & numberDone);

private:
  /** Default Constructor, forbiden for some reason.*/
  BonChooseVariable ();

  /// CbcModel, used to get status of search
  CbcModel* cbc_model_;

  /** Flag indicating whether we don't want to mix strong branching
   *  and pseudo costs during the decision which variable to branch
   *  on */
  bool only_pseudo_when_trusted_;

  /** Number of variables put into the list because there were not
   *  trusted */
  int number_not_trusted_;

  // ToDo: Make this options
  /** @name Algoirithmic options */
  //@{
  /** maxmin weight in branching decision when no solution has been
   *  found yet */
  double maxmin_crit_no_sol_;
  /** maxmin weight in branching decision when no solution has been
   *  found yet */
  double maxmin_crit_have_sol_;
  /** fraction of branching candidates that are not trusted yet */
  double setup_pseudo_frac_;
  /** number of times a branch has to happen so that it is trusted in
   *  setupList */
  int numberBeforeTrustedList_;
  /** number of strong branching points at root node */
  int numberStrongRoot_;
  /** backup of numberStrong_ before Root node solve */
  int numberStrongBackup_;
  //@}

  double maxminCrit() const;

  /** detecting if this is root node */
  bool isRootNode() const;

  GuessHeuristic* guessHeuristic_;
};

}
#endif
