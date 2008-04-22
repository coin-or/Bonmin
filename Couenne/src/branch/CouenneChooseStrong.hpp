/*
 * Name:    CouenneChooseStrong.hpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching object for Couenne
 *
 * (C) Carnegie-Mellon University, 2006.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECHOOSESTRONG_HPP
#define COUENNECHOOSESTRONG_HPP

#include "BonChooseVariable.hpp"

// Forward declaration
class CouenneProblem;

namespace Bonmin {

  class CouenneChooseStrong : public BonChooseVariable  {

  public:

    enum pseudocostMult {infeasibility, interval, projectDist};

    /// Constructor from solver (so we can set up arrays etc)
    CouenneChooseStrong (BabSetupBase& b, CouenneProblem* problem);

    /// Copy constructor
    CouenneChooseStrong (const CouenneChooseStrong &);

    /// Assignment operator
    CouenneChooseStrong & operator= (const CouenneChooseStrong& rhs);

    /// Clone
    virtual OsiChooseVariable * clone() const;

    /// Destructor
    virtual ~CouenneChooseStrong ();

    // Sets up strong list and clears all if initialize is true.
    // Returns number of infeasibilities.
    virtual int setupList (OsiBranchingInformation *info, bool initialize);

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
    virtual int doStrongBranching( OsiSolverInterface * solver, 
				   OsiBranchingInformation *info,
				   int numberToDo, int returnCriterion);

    /// Add list of options to be read from file
    static void registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions);

  private:

    /** Default Constructor, forbidden for some reason.*/
    CouenneChooseStrong ();
    
    /// Pointer to the associated MINLP problem
    CouenneProblem *problem_;

    /// multiplier type for pseudocost
    enum pseudocostMult pseudoMultType_;

    /// should we update the pseudocost multiplier with the distance
    /// between the LP point and the solution of the resulting
    /// branches' LPs?
    bool pseudoUpdateLP_;
  };
}
#endif
