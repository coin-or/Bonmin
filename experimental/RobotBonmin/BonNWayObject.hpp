// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Pierre Bonami 13/10/2011 copy and modify CbcNWay to adapt to use an OsiObject for Bonmin
//
// Edwin 11/9/2009-- carved out of CbcBranchActual

/** Define an n-way class for variables.
    Only valid value is one at UB others at LB
    Normally 0-1
*/
#ifndef BonNWayObject_H
#define BonNWayObject_H
#include "OsiBranchingObject.hpp"
#include "CbcConsequence.hpp"
#include <list>

namespace Bonmin {
class n_way_consequences {
public:
  std::vector<int> indices;

   n_way_consequences(): indices(){
   }

   n_way_consequences(const n_way_consequences& rhs): 
        indices(rhs.indices){
   }

   n_way_consequences& operator=(const n_way_consequences& rhs){
     if(this != &rhs){
       indices = rhs.indices;
     }
     return *this;
   }

  n_way_consequences * clone() const{
     return new n_way_consequences(*this);
  }

  void applyToSolver(OsiSolverInterface * solver, int state) const{
     if(state < 0) return;
     int num_one = 0;
     for(size_t i = 0 ; i < indices.size() ; i++){
         if(solver->getColLower()[indices[i]] > 0.9) num_one++;
         solver->setColUpper(indices[i],solver->getColLower()[indices[i]]);
     }
     assert(num_one == 0);
  }

};

class BonNWayObject : public OsiObject {

public:

    // Default Constructor
    BonNWayObject ();

    /** Useful constructor (which are matrix indices)
    */
    BonNWayObject (int numberMembers,
             const int * which, int identifier);

    // Copy constructor
    BonNWayObject ( const BonNWayObject &);

    /// Clone
    virtual OsiObject * clone() const;

    /// Assignment operator
    BonNWayObject & operator=( const BonNWayObject& rhs);

    /// Destructor
    virtual ~BonNWayObject ();

    /// Set up a consequence for a single member
    void setConsequence(int iMember, const n_way_consequences & consequence);

    /// Applies a consequence for a single member
    void applyConsequence(OsiSolverInterface * solver,
                          int iSequence, int state) const;

    /// Infeasibility - large is 0.5 (and 0.5 will give this)
    virtual double infeasibility(const OsiBranchingInformation * info,
                                 int &preferredWay) const;

    /// This looks at solution and sets bounds to contain solution
    virtual double feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const;

    /// Creates a branching object
    virtual OsiBranchingObject * createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const;

    /// Number of members
    inline size_t numberMembers() const {
        return members_.size();
    }

    /// Members (indices in range 0 ... numberColumns-1)
    inline const int * members() const {
        return members_.data();
    }

    void set_bounds(std::vector<double> & bounds) const{
      bounds_ = bounds;
    }

    void make_quick(){quicky_ = true;}

    void set_only_frac_branches(int depth){
      only_frac_branch_ = depth;
    }
private:
    /// data

    /// Members (indices in range 0 ... numberColumns-1)
    std::vector<int> members_;
    /// Consequences (normally NULL)
    n_way_consequences ** consequence_;

    /// Bounds on the members
    mutable std::vector<double> bounds_;

    /// Quicky only branch up on variables with non zero value
    bool quicky_;
    /// Only branch on fractional variables (last branch puts all of them to 0)
    int only_frac_branch_;
};
/** N way branching Object class.
    Variable is number of set.
 */
class BonNWayBranchingObject : public OsiBranchingObject {

public:

    // Default Constructor
    BonNWayBranchingObject ();

    /** Useful constructor - order had matrix indices
        way_ -1 corresponds to setting first, +1 to second, +3 etc.
        this is so -1 and +1 have similarity to normal
    */
    BonNWayBranchingObject (OsiSolverInterface * solver,  const BonNWayObject * nway,
                            const std::vector<int>& order, const std::list<int>& skipped);

    // Copy constructor
    BonNWayBranchingObject ( const BonNWayBranchingObject &);

    // Assignment operator
    BonNWayBranchingObject & operator=( const BonNWayBranchingObject& rhs);

    /// Clone
    virtual OsiBranchingObject * clone() const;

    // Destructor
    virtual ~BonNWayBranchingObject ();

    using OsiBranchingObject::branch ;
    /// Does next branch and updates state
    virtual double branch(OsiSolverInterface * solver);

    /// Branch on specified  state
    virtual double branch(OsiSolverInterface * solver, int state);

    /** The number of branch arms created for this branching object
    */
    virtual int numberBranches() const {
        return static_cast<int>(order_.size()) + (!skipped_.empty());
    }
    /// Is this a two way object (-1 down, +1 up)
    virtual bool twoWay() const {
        return false;
    }

    inline int var_branched_on(){
      if(branchIndex_ > 0)
        return object_->members()[order_[branchIndex_ - 1]];
      return -1;
    }

    inline int seq_branched_on(){
      if(branchIndex_ > 0)
      return order_[branchIndex_ - 1];
      return -1;
    }

    inline int state_branched_on(){
      return branchIndex_ - 1;
    }
private:
    /// Points back to object
    const BonNWayObject * object_;
    /// order of branching 
    std::vector<int> order_;
    /// Is only branching on a subset of variables (has to do a last branch with all variables in order set to 0)
    std::list<int> skipped_;
};

}

#endif
