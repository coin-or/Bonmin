// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <climits>


#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "BonNWayObject.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//#define FULL_PRINT
namespace Bonmin {
//##############################################################################

// Default Constructor
BonNWayObject::BonNWayObject ()
        : OsiObject(),
        members_(),
        consequence_(NULL),
        quicky_(false),
        only_frac_branch_(0)
{
}

// Useful constructor (which are integer indices)
BonNWayObject::BonNWayObject (int numberMembers,
                  const int * which, int identifier):
       OsiObject(),
       members_(numberMembers),
        quicky_(false),
        only_frac_branch_(0)//INT_MAX)
{
    if (numberMembers) {
        memcpy(members_.data(), which, numberMembers*sizeof(int));
    } else {
    }
    consequence_ = NULL;
}

// Copy constructor
BonNWayObject::BonNWayObject ( const BonNWayObject & rhs)
        : OsiObject(rhs), members_(rhs.members_), quicky_(rhs.quicky_), only_frac_branch_(rhs.only_frac_branch_)
{
    consequence_ = NULL;
    if (members_.size()) {
        if (rhs.consequence_) {
            consequence_ = new n_way_consequences* [members_.size()];
            for (size_t i = 0; i < members_.size(); i++) {
                if (rhs.consequence_[i])
                    consequence_[i] = rhs.consequence_[i]->clone();
                else
                    consequence_[i] = NULL;
            }
        }
    } else {
    }
}

// Clone
OsiObject *
BonNWayObject::clone() const
{
    return new BonNWayObject(*this);
}

// Assignment operator
BonNWayObject &
BonNWayObject::operator=( const BonNWayObject & rhs)
{
    if (this != &rhs) {
        OsiObject::operator=(rhs);
        members_ = rhs.members_;
        quicky_ = rhs.quicky_;
        only_frac_branch_ = rhs.only_frac_branch_;
        if (consequence_) {
            for (size_t i = 0; i < members_.size(); i++)
                delete consequence_[i];
            delete [] consequence_;
            consequence_ = NULL;
        }
        if (rhs.consequence_) {
            consequence_ = new n_way_consequences* [members_.size()];
            for (size_t i = 0; i < members_.size(); i++) {
                if (rhs.consequence_[i])
                    consequence_[i] = rhs.consequence_[i]->clone();
                else
                    consequence_[i] = NULL;
            }
        }
    }
    return *this;
}

// Destructor
BonNWayObject::~BonNWayObject ()
{
    if (consequence_) {
        for (size_t i = 0; i < members_.size(); i++)
            delete consequence_[i];
        delete [] consequence_;
    }
}
// Set up a consequence for a single member
void
BonNWayObject::setConsequence(int iMember, const n_way_consequences& consequence)
{
    if (!consequence_) {
        consequence_ = new n_way_consequences* [members_.size()];
        for (size_t i = 0; i < members_.size(); i++)
            consequence_[i] = NULL;
    }
    consequence_[iMember] = consequence.clone();
}

// Applies a consequence for a single member
void
BonNWayObject::applyConsequence(OsiSolverInterface * solver, int iSequence, int state) const
{
    assert (state == -9999 || state == 9999);
    if (consequence_) {
        n_way_consequences* consequence = consequence_[iSequence];
        if (consequence)
            consequence->applyToSolver(solver, state);
    }
}
double
BonNWayObject::infeasibility(const OsiBranchingInformation * info,
                       int &preferredWay) const
{
    int numberUnsatis = 0;
    size_t j;
    const double * solution = info->solution_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double largestValue = 0.0;

    double integerTolerance = info->integerTolerance_;

    for (j = 0; j < members_.size(); j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        double distance = CoinMin(value - lower[iColumn], upper[iColumn] - value);
        if (distance > integerTolerance) {
            numberUnsatis++;
            largestValue = CoinMax(distance, largestValue);
        }
    }
    preferredWay = 1;
    if (numberUnsatis) {
        return largestValue;
    } else {
        return 0.0; // satisfied
    }
}

// This looks at solution and sets bounds to contain solution
double
BonNWayObject::feasibleRegion(OsiSolverInterface * solver,
                        const OsiBranchingInformation * info) const
{
    size_t j;
    const double * solution = info->solution_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double r_val = 0;
    for (j = 0; j < members_.size(); j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (value >= upper[iColumn] - integerTolerance) {
            r_val += value - upper[iColumn];
            solver->setColLower(iColumn, upper[iColumn]);
        } else {
            assert (value <= lower[iColumn] + integerTolerance);
            r_val += value - lower[iColumn];
            solver->setColUpper(iColumn, lower[iColumn]);
        }
    }
    return r_val;
}

OsiBranchingObject *
BonNWayObject::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{

    const double * solution = info->solution_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double cutoff = info->cutoff_;
    std::vector<int> list;
    list.reserve(members_.size());
    std::vector<double> sort;
    sort.reserve(members_.size());

    int n_skipped = 0;
    std::list<int> skipped;

    for (size_t j = 0; j < members_.size(); j++) {
        const int& iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (upper[iColumn] > lower[iColumn]) {
        if(bounds_.size() && bounds_[j] > cutoff) {
          //printf("Fathom branch on bound %i : %g > %g\n", j, bounds_[j], cutoff);
          continue;
        }
        if( (quicky_ || info->depth_ > only_frac_branch_ ) && fabs(solution[iColumn] - lower[iColumn]) < integerTolerance){
            if(info->depth_ > only_frac_branch_) {
               n_skipped++;
               //printf("Skipping variable %i\n", iColumn);
               skipped.push_back(static_cast<int>(j));
               continue;
            }
          }
            double distance = upper[iColumn] - value;
            list.push_back(static_cast<int>(j));
            sort.push_back(distance);
        }
    }

    if(n_skipped == 1) {//False subset put back the missing since branching up allows applying consequences
           const int& iColumn = members_[skipped.front()];
           double value = solution[iColumn];
           value = CoinMax(value, lower[iColumn]);
           value = CoinMin(value, upper[iColumn]);
           double distance = upper[iColumn] - value;
           list.push_back(skipped.front());
           sort.push_back(distance);
           skipped.pop_front();
           n_skipped --;
    }

    // sort
    CoinSort_2(sort.data(), sort.data() + sort.size(), list.data());
    // create object
    OsiBranchingObject * branch;

    //if(n_skipped) printf("Creating branch n_skipped is %i numberFree %i\n", n_skipped, numberFree);
    branch = new BonNWayBranchingObject(solver, this, list, skipped);
    return branch;
}


// Default Constructor
BonNWayBranchingObject::BonNWayBranchingObject()
        : OsiBranchingObject(), object_(NULL),
          order_(), skipped_()
{
}

// Useful constructor
BonNWayBranchingObject::BonNWayBranchingObject (OsiSolverInterface * solver,
        const BonNWayObject * nway,
        const std::vector<int>& order, const std::list<int> &skipped)
        : OsiBranchingObject(solver, 0.5), order_(order), skipped_(skipped)
{
    numberBranches_ = static_cast<int>(order_.size()) + (!skipped.empty());
    object_ = nway;
}

// Copy constructor
BonNWayBranchingObject::BonNWayBranchingObject ( const BonNWayBranchingObject & rhs) :
           OsiBranchingObject(rhs), object_(rhs.object_), 
           order_(rhs.order_), skipped_(rhs.skipped_)
{
    object_ = rhs.object_;
}

// Assignment operator
BonNWayBranchingObject &
BonNWayBranchingObject::operator=( const BonNWayBranchingObject & rhs)
{
    if (this != &rhs) {
        OsiBranchingObject::operator=(rhs);
        object_ = rhs.object_;
        order_ = rhs.order_;
        skipped_ = rhs.skipped_;
    }
    return *this;
}
OsiBranchingObject *
BonNWayBranchingObject::clone() const
{
    return (new BonNWayBranchingObject(*this));
}


// Destructor
BonNWayBranchingObject::~BonNWayBranchingObject ()
{
}

double
BonNWayBranchingObject::branch(OsiSolverInterface * solver)
{
    int which = branchIndex_;
    branchIndex_++;
    if(!skipped_.empty() && branchIndex_ == static_cast<int>(order_.size())){//We are branching on a subset last branch fixes all in subset to 0
      which = -1; // Simply done by setting which to a dummy value;
    }
    return branch(solver, which);
}

double
BonNWayBranchingObject::branch(OsiSolverInterface * solver, int which)
{
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const int * members = object_->members();
    for (size_t j = 0; j < order_.size(); j++) {
        const int& iSequence = order_[j];
        const int& iColumn = members[iSequence];
        if (j != which) {
            solver->setColUpper(iColumn, lower[iColumn]);
            assert (lower[iColumn] > -1.0e20);
        } else {
            solver->setColLower(iColumn, upper[iColumn]);
#ifdef FULL_PRINT
            printf("Up Fix %d to %g\n", iColumn, upper[iColumn]);
#endif
            assert (upper[iColumn] < 1.0e20);
            // apply any consequences
            object_->applyConsequence(solver, iSequence, 9999);
        }
    }
    if(which != -1){
       for(std::list<int>::iterator k = skipped_.begin() ; k != skipped_.end() ; k++){
           assert(upper[members[*k]] > lower[members[*k]] + 0.5);
           solver->setColUpper(members[*k], lower[members[*k]]);
       }
    }
    return 0.0;
}

}

