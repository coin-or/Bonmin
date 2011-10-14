// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/9/2009-- carved out of CbcBranchActual

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

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
        consequence_(NULL)
{
}

// Useful constructor (which are integer indices)
BonNWayObject::BonNWayObject (int numberMembers,
                  const int * which, int identifier):
       OsiObject(),
       members_(numberMembers)
{
    if (numberMembers) {
        memcpy(members_.data(), which, numberMembers*sizeof(int));
    } else {
    }
    consequence_ = NULL;
}

// Copy constructor
BonNWayObject::BonNWayObject ( const BonNWayObject & rhs)
        : OsiObject(rhs), members_(rhs.members_)
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
    int numberFree = 0;
    size_t j;

    const double * solution = info->solution_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double cutoff = info->cutoff_;
    int * list = new int[members_.size()];
    double * sort = new double[members_.size()];

    for (j = 0; j < members_.size(); j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (upper[iColumn] > lower[iColumn]) {
        if(bounds_.size() && bounds_[j] > cutoff) {
          printf("Fathom branch on bound %i : %g > %g\n", j, bounds_[j], cutoff);
          continue;
        }
            double distance = upper[iColumn] - value;
            list[numberFree] = static_cast<int>(j);
            sort[numberFree++] = distance;
        }
    }
    if(numberFree ==0){//Could be fathomed but impossible here.... create one branch
      for (j = 0; j < members_.size(); j++) {
        int iColumn = members_[j];
        double value = solution[iColumn];
        value = CoinMax(value, lower[iColumn]);
        value = CoinMin(value, upper[iColumn]);
        if (upper[iColumn] > lower[iColumn]) {
            double distance = upper[iColumn] - value;
            list[numberFree] = static_cast<int>(j);
            sort[numberFree++] = distance;
            break;
        }
    }
    }
    assert (numberFree);
    // sort
    CoinSort_2(sort, sort + numberFree, list);
    // create object
    OsiBranchingObject * branch;
    branch = numberFree ? new BonNWayBranchingObject(solver, this, numberFree, list) : NULL;
    delete [] list;
    delete [] sort;
    return branch;
}

// Default Constructor
BonNWayBranchingObject::BonNWayBranchingObject()
        : OsiBranchingObject()
{
    order_ = NULL;
    object_ = NULL;
    numberInSet_ = 0;
    way_ = 0;
}

// Useful constructor
BonNWayBranchingObject::BonNWayBranchingObject (OsiSolverInterface * solver,
        const BonNWayObject * nway,
        int number, const int * order)
        : OsiBranchingObject(solver, 0.5)
{
    numberBranches_ = number;
    order_ = new int [number];
    object_ = nway;
    numberInSet_ = number;
    memcpy(order_, order, number*sizeof(int));
    way_ = -1;
}

// Copy constructor
BonNWayBranchingObject::BonNWayBranchingObject ( const BonNWayBranchingObject & rhs) : OsiBranchingObject(rhs)
{
    numberInSet_ = rhs.numberInSet_;
    object_ = rhs.object_;
    way_ = rhs.way_;
    if (numberInSet_) {
        order_ = new int [numberInSet_];
        memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
    } else {
        order_ = NULL;
    }
}

// Assignment operator
BonNWayBranchingObject &
BonNWayBranchingObject::operator=( const BonNWayBranchingObject & rhs)
{
    if (this != &rhs) {
        OsiBranchingObject::operator=(rhs);
        object_ = rhs.object_;
        way_ = rhs.way_;
        delete [] order_;
        numberInSet_ = rhs.numberInSet_;
        if (numberInSet_) {
            order_ = new int [numberInSet_];
            memcpy(order_, rhs.order_, numberInSet_*sizeof(int));
        } else {
            order_ = NULL;
        }
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
    delete [] order_;
}

double
BonNWayBranchingObject::branch(OsiSolverInterface * solver)
{
    int which = branchIndex_;
    branchIndex_++;
    return branch(solver, which);
}

double
BonNWayBranchingObject::branch(OsiSolverInterface * solver, int which)
{
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const int * members = object_->members();
    for (int j = 0; j < numberInSet_; j++) {
        int iSequence = order_[j];
        int iColumn = members[iSequence];
        if (j != which) {
            solver->setColUpper(iColumn, lower[iColumn]);
            //model_->solver()->setColLower(iColumn,lower[iColumn]);
            assert (lower[iColumn] > -1.0e20);
            // apply any consequences
            object_->applyConsequence(solver, iSequence, -9999);
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
    return 0.0;
}

}

