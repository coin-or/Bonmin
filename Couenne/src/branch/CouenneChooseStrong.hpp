/*
 * Name:    CouenneChooseStrong.hpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching object for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECHOOSESTRONG_HPP
#define COUENNECHOOSESTRONG_HPP

#include "BonChooseVariable.hpp"
#include "CouenneJournalist.hpp"


// Forward declaration
class CouenneProblem;

namespace Bonmin {

  class CouenneChooseStrong : public BonChooseVariable  {

  public:

    /// Constructor from solver (so we can set up arrays etc)
    CouenneChooseStrong (BabSetupBase& b, CouenneProblem* problem, JnlstPtr jnlst);

    /// Copy constructor
    CouenneChooseStrong (const CouenneChooseStrong &);

    /// Assignment operator
    CouenneChooseStrong & operator= (const CouenneChooseStrong& rhs);

    /// Clone
    virtual OsiChooseVariable * clone() const;

    /// Destructor
    virtual ~CouenneChooseStrong ();

    /// Sets up strong list and clears all if initialize is true.
    /// Returns number of infeasibilities.
    virtual int setupList (OsiBranchingInformation *info, bool initialize);

    /**  This is a utility function which does strong branching on a
	 list of objects and stores the results in OsiHotInfo.objects.
	 On entry the object sequence is stored in the OsiHotInfo
	 object and maybe more.
	 It returns -
	 -1 - one branch was infeasible both ways
	 0 - all inspected - nothing can be fixed
	 1 - all inspected - some can be fixed                         (returnCriterion==0)
	 2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
	 3 - returning because max time
    */
    virtual int doStrongBranching( OsiSolverInterface * solver, 
				   OsiBranchingInformation *info,
				   int numberToDo, int returnCriterion);

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

  private:

    /** Default Constructor, forbidden for some reason.*/
    CouenneChooseStrong ();
    
  protected:

    /// Pointer to the associated MINLP problem
    CouenneProblem *problem_;

    /// should we update the pseudocost multiplier with the distance
    /// between the LP point and the solution of the resulting
    /// branches' LPs? If so, this only happens in strong branching
    bool pseudoUpdateLP_;

    /// pointer to journalist for detailed information
    JnlstPtr jnlst_;
  };
}

#endif
