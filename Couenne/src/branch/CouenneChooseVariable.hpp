/*
 * Name:    CouenneChooseVariable.hpp
 * Author:  Pietro Belotti
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECHOOSEVARIABLE_HPP
#define COUENNECHOOSEVARIABLE_HPP

#include <OsiChooseVariable.hpp>


/** \brief Choose a variable for branching
 */

class CouenneChooseVariable: public OsiChooseVariable {

public:

  /// Default Constructor 
  CouenneChooseVariable ();

  /// Constructor from solver (so we can set up arrays etc)
  CouenneChooseVariable (const OsiSolverInterface *solver);

  /// Copy constructor 
  CouenneChooseVariable (const CouenneChooseVariable &);

  /// Assignment operator 
  CouenneChooseVariable &operator= (const CouenneChooseVariable& rhs);

  /// Clone
  virtual OsiChooseVariable *clone() const;

  /// Destructor 
  virtual ~CouenneChooseVariable ();

  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. 
      If returns -1 then has worked out node is infeasible!
  */

  virtual int setupList (OsiBranchingInformation *info, bool initialize);


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

  virtual int chooseVariable (OsiSolverInterface * solver, 
			      OsiBranchingInformation *info, 
			      bool fixVariables);


  /// Given a candidate fill in useful information e.g. estimates
  virtual void updateInformation( const OsiBranchingInformation *info,
				  int branch, OsiHotInfo * hotInfo);
};

#endif
