// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 11/21/2006

#include "BonFilterWarmStart.hpp"


namespace Bonmin
{

  FilterWarmStart:: FilterWarmStart(const fint xSize /*= 0*/,
      const real* xArray /*= NULL*/,
      const fint lamSize /*= 0*/,
      const real* lamArray /*= NULL*/,
      const fint lwsSize /*= 0*/,
      const fint *lwsArray /*= NULL*/,
      const fint istat[14] /*= def_istat*/)
      :
      CoinWarmStartBasis(),
      xSize_(xSize),
      xArray_(NULL),
      tempxArray_(NULL),
      lamSize_(lamSize),
      lamArray_(NULL),
      templamArray_(NULL),
      lwsSize_(lwsSize),
      lwsArray_(NULL)
  {
    DBG_ASSERT(xSize > 0 || !xArray);
    if (xSize_ > 0) {
      xArray_ = new real[xSize];
      DBG_ASSERT(xArray);
      CoinCopyN(xArray, xSize, xArray_);
    }
    DBG_ASSERT(lamSize > 0 || !lamArray);
    if (lamSize_ > 0) {
      lamArray_ = new real[lamSize];
      DBG_ASSERT(lamArray);
      CoinCopyN(lamArray, lamSize, lamArray_);
    }
    DBG_ASSERT(lwsSize > 0 || !lwsArray);
    if (lwsSize_ > 0) {
      lwsArray_ = new fint[lwsSize];
      DBG_ASSERT(lwsArray);
      CoinCopyN(lwsArray, lwsSize, lwsArray_);
    }
    for (int i = 0 ; i < 14 ; i ++) {
      istat_[i] = istat[i];
    }
  }

  /* Copy constructor */
  FilterWarmStart::FilterWarmStart(const FilterWarmStart & other)
      :
      CoinWarmStartBasis(other),
      xSize_(other.xSize_),
      xArray_(NULL),
      tempxArray_(other.tempxArray_),
      lamSize_(other.lamSize_),
      lamArray_(NULL),
      templamArray_(other.templamArray_),
      lwsSize_(other.lwsSize_),
      lwsArray_(NULL)
  {
    DBG_ASSERT(other.xSize_ > 0 || !other.xArray_);
    if (xSize_ > 0) {
      xArray_ = new real[xSize_];
      DBG_ASSERT(other.xArray_);
      CoinCopyN(other.xArray_, xSize_, xArray_);
    }
    DBG_ASSERT(lamSize_ > 0 || !lamArray_);
    if (lamSize_ > 0) {
      lamArray_ = new real[lamSize_];
      DBG_ASSERT(other.lamArray_);
      CoinCopyN(other.lamArray_, lamSize_, lamArray_);
    }
    DBG_ASSERT(lwsSize_ > 0 || !lwsArray_);
    if (lwsSize_ > 0) {
      lwsArray_ = new fint[lwsSize_];
      DBG_ASSERT(other.lwsArray_);
      CoinCopyN(other.lwsArray_, lwsSize_, lwsArray_);
    }
    for (int i = 0 ; i < 14 ; i ++) {
      istat_[i] = other.istat_[i];
    }
  }

  FilterWarmStart::~FilterWarmStart()
  {
    delete [] xArray_;
    delete [] lamArray_;
    delete [] lwsArray_;
  }

  CoinWarmStartDiff *
  FilterWarmStart::generateDiff(const CoinWarmStart * const oldOne) const
  {
    const FilterWarmStart * old =
      dynamic_cast<const FilterWarmStart*> (oldOne);

    if (xSize_ != old->xSize_ || lamSize_ != old->lamSize_
        || lwsSize_ != old->lwsSize_) {
      throw CoinError("Can not make difference for warm starts of differnet sizes",
          "generateDiffs",
          "FilterWarmStart");
    }

    FilterWarmStartDiff* diff =
      new FilterWarmStartDiff(xSize_, xArray_, lamSize_, lamArray_, lwsSize_);

    for (fint i = 0 ; i < lwsSize_ ; i++) {
      if (lwsArray_[i] != old->lwsArray_[i]) {
        diff->differences.push_back(FilterWarmStartDiff::OneDiff(i, lwsArray_[i] - old->lwsArray_[i]));
      }
    }

    diff->differences.resize(diff->differences.size());

    for (int i = 0 ; i < 14 ; i++) {
      diff->istat_[i] = istat_[i];
    }
    return diff;
  }


  void
  FilterWarmStart::applyDiff(const CoinWarmStartDiff * diff)
  {

    const FilterWarmStartDiff * diffF =
      dynamic_cast<const FilterWarmStartDiff  *>(diff);
    DBG_ASSERT(diffF != NULL);

    tempxArray_ = diffF->xArray_;
    templamArray_ = diffF->lamArray_;

    fint end = diffF->differences.size();
    for (fint i = 0 ; i < end ; i++) {
      lwsArray_[diffF->differences[i].first] += diffF->differences[i].second;
    }

    for (int i = 0 ; i < 14 ; i++)
      istat_[i] = diffF->istat_[i];
  }

  FilterWarmStartDiff::FilterWarmStartDiff(fint xSize,
      real* xArray,
      fint lamSize,
      real* lamArray,
      fint capacity)
      :
      CoinWarmStartBasisDiff()
  {
    xSize_ = xSize;
    xArray_ = new real[xSize];
    CoinCopyN(xArray, xSize, xArray_);
    lamSize_ = lamSize;
    lamArray_ = new real[lamSize];
    CoinCopyN(lamArray, lamSize, lamArray_);

    differences.reserve(capacity);
  }

  FilterWarmStartDiff::~FilterWarmStartDiff()
  {
    delete [] xArray_;
    delete [] lamArray_;
  }

  void
  FilterWarmStart::flushPoint()
  {
    delete [] xArray_;
    delete [] lamArray_;
    delete [] lwsArray_;

    xArray_ = NULL;
    lamArray_ = NULL;
    lwsArray_ = NULL;
  }

  void
  FilterWarmStartDiff::flushPoint()
  {
    delete [] xArray_;
    delete [] lamArray_;

    xArray_ = NULL;
    lamArray_ = NULL;
  }

  FilterSolver::fint
  FilterWarmStart::def_istat[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


} /* End namespace Bonmin */
