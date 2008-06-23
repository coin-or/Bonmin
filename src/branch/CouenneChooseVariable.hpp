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
#include "BonBabInfos.hpp"
#include "CouenneJournalist.hpp"

class CouenneProblem;

/** \brief Choose a variable for branching
 */

class CouenneChooseVariable: public OsiChooseVariable {

public:

  /// Default Constructor
  CouenneChooseVariable ();

  /// Constructor from solver (so we can set up arrays etc)
  CouenneChooseVariable (const OsiSolverInterface *, CouenneProblem *, JnlstPtr jnlst);

  /// Copy constructor 
  CouenneChooseVariable (const CouenneChooseVariable &);

  /// Assignment operator 
  CouenneChooseVariable &operator= (const CouenneChooseVariable &);

  /// Clone
  virtual OsiChooseVariable *clone() const
  {return new CouenneChooseVariable (*this);}

  /// Destructor 
  virtual ~CouenneChooseVariable () {}

  /// Sets up strong list and clears all if initialize is true.
  /// Returns number of infeasibilities.  If returns -1 then has
  /// worked out node is infeasible!
  virtual int setupList (OsiBranchingInformation *, bool);

  /// Returns true if solution looks feasible against given objects
  virtual bool feasibleSolution (const OsiBranchingInformation * info,
				 const double * solution,
				 int numberObjects,
				 const OsiObject ** objects);

  /// choose object to branch based on earlier setup
  virtual int chooseVariable (OsiSolverInterface * solver, 
			      OsiBranchingInformation *info, 
			      bool fixVariables);

  /// Add list of options to be read from file
  static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

protected:

  /// Pointer to the associated MINLP problem
  CouenneProblem *problem_;

  /// journalist for detailed debug information
  JnlstPtr jnlst_;
};

#endif
