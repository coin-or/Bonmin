/*
 * Name:    CouenneChooseVariable.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECHOOSEVARIABLE_HPP
#define COUENNECHOOSEVARIABLE_HPP

#include "OsiChooseVariable.hpp"
#include "CouenneProblem.hpp"
#include "BonBabSetupBase.hpp"

// Forward declaration
class CbcModel;

/** \brief Choose a variable for branching
 */

class CouenneChooseVariable: public OsiChooseVariable {

public:

  /// Default Constructor 
  CouenneChooseVariable ();

  /// Constructor from solver (so we can set up arrays etc)
  CouenneChooseVariable (const OsiSolverInterface *, CouenneProblem *);

  /// Copy constructor 
  CouenneChooseVariable (const CouenneChooseVariable &);

  /// Assignment operator 
  CouenneChooseVariable &operator= (const CouenneChooseVariable &);

  /// Clone
  virtual OsiChooseVariable *clone() const;

  /// Destructor 
  virtual ~CouenneChooseVariable ();

  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. 
      If returns -1 then has worked out node is infeasible!
  */

  virtual int setupList (OsiBranchingInformation *, bool);

  /** Choose a variable
      Returns:
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from bestObjectIndex() and bestWhichWay()
     We can pick up a forced branch (can change bound) from firstForcedObjectIndex() 
     and firstForcedWhichWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
     If fixVariables is true then 2,3,4 are all really same as problem changed
  */

  virtual int chooseVariable (OsiSolverInterface *, 
			      OsiBranchingInformation *, 
			      bool);

  /// Return pointer to associated MINLP problem
  CouenneProblem *Problem () const 
  {return problem_;}

protected:

  /// Pointer to the associated MINLP problem
  CouenneProblem *problem_;
};

class CouenneChooseStrong : public OsiChooseStrong  {

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
  CouenneChooseStrong (Bonmin::BabSetupBase& b, CouenneProblem* si);

  /// Copy constructor 
  CouenneChooseStrong (const CouenneChooseStrong &);
   
  /// Assignment operator 
  CouenneChooseStrong & operator= (const CouenneChooseStrong& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~CouenneChooseStrong ();

  /** Helper functions for setupList and chooseVariable */
  double maxminCrit(const OsiBranchingInformation* info) const;
  void computeMultipliers(double& upMult, double& downMult) const;
  double computeUsefulness(const double MAXMIN_CRITERION,
			   const double upMult, const double dowMult,
			   const double value,
			   const OsiObject* object, int i,
			   double& value2) const;

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
#if 0
  /// Given a branch fill in useful information e.g. estimates
  virtual void updateInformation( int whichObject, int branch, 
				  double changeInObjective, double changeInValue,
				  int status);
#endif

  /** Method for setting CbcModel, which is used to get statusOfSearch */
  void setCbcModel(CbcModel* cbc_model)
  {
    cbc_model_ = cbc_model;
  }

  void setOnlyPseudoWhenTrusted(bool only_pseudo_when_trusted)
  {
    only_pseudo_when_trusted_ = only_pseudo_when_trusted;
  }

protected:

  /// Holding on the a pointer to the journalist
  SmartPtr<Journalist> jnlst_;

  /// verbosity level
  int bb_log_level_;

private:
  /** Default Constructor, forbiden for some reason.*/
  CouenneChooseStrong ();

  /// CbcModel, used to get status of search
  CbcModel* cbc_model_;

  /** Flag indicating whether we don't want to mix strong branching
   *  and pseudo costs during the decision which variable to branch
   *  on */
  bool only_pseudo_when_trusted_;

  /** Number of variables put into the list because there were not
   *  trusted */
  int number_not_trusted_;

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

  /** detecting if this is root node */
  bool isRootNode(const OsiBranchingInformation *info) const;

  /// Pointer to the associated MINLP problem
  CouenneProblem *problem_;
};


#endif
