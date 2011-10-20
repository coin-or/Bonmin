// Copyright (C) 2006, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonNWayChoose_H
#define BonNWayChoose_H

#include "OsiChooseVariable.hpp"
#include "BonBabSetupBase.hpp"
#include "BonNWayObject.hpp"

#define OLD_USEFULLNESS

namespace Bonmin
{

  /** This class chooses a variable to branch on

      This is the base class for the branching rules in Bonmin (inherits
      from OsiChooseVariable). This class implements a simple strong branching algorithm where the changes in the objective
      value induced by branching on a specific object are estimated with the pure virtual function fill_changes.
  */

  class BonNWayChoose : public OsiChooseVariable
  {
  protected:

  public:

    /// Constructor from solver (so we can set up arrays etc)
    BonNWayChoose (BabSetupBase& b, const OsiSolverInterface* solver);

    /// Copy constructor
    BonNWayChoose (const BonNWayChoose &);

    /// Assignment operator
    BonNWayChoose & operator= (const BonNWayChoose& rhs);

    /// Clone
    virtual OsiChooseVariable * clone() const;

    /// Destructor
    virtual ~BonNWayChoose ();

    /** Helper function for setupList and chooseVariable, compute usefullness of an nway object */
    double compute_usefulness(int objectIndex, const OsiBranchingInformation * info) const;

    double compute_usefulness(const OsiBranchingInformation * info,
                size_t n, const int * vars, const std::vector<double> &bounds, const std::vector<double> &unit_changes) const;
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
  virtual int doStrongBranching( OsiSolverInterface * solver, 
				 OsiBranchingInformation *info,
				 int iObject, double * saveLow, double * saveUp, double &score);

 
  /** Register class options.*/
  static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

  private:
    /** Default Constructor, forbiden for some reason.*/
    BonNWayChoose ();

    /** depth of strong-branching.*/
    int br_depth_;
    /** Do we fix?*.*/
    Ipopt::Index do_fixings_;
    /** Cutoff multiplier for infeasibles*/
    double cutoff_multiplier_;
    /** Trust value for the pseudo costs.*/
    double pseudocost_trust_value_;
    /** Global time limit for algorithm. */
    double time_limit_;

    /** Starting time of algorithm.*/
    double start_time_;

    /// Start of nway objects in array
    int start_nway_;

    /// log level
    int log_;
    typedef std::vector< std::vector<double> > full_mat;

    full_mat bounds_;

    full_mat unit_changes_;

    std::vector < std::vector<int> > num_ps_costs_;

    std::vector<int> num_eval_;

    int geo_means_;
  };

}
#endif
