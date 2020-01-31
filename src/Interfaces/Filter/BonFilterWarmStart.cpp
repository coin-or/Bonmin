// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
    CoinWarmStartPrimalDual(xSize, lamSize, xArray, lamArray),
    CoinWarmStartBasis(),
    lwsSize_(lwsSize),
    lwsArray_(NULL),
    empty_(false)
  {
    assert(lwsSize > 0 || !lwsArray);
    if (lwsSize_ > 0) {
      lwsArray_ = new fint[lwsSize];
      assert(lwsArray);
      CoinCopyN(lwsArray, lwsSize, lwsArray_);
    }
    for (int i = 0 ; i < 14 ; i ++) {
      istat_[i] = istat[i];
    }
  }

  /* Copy constructor */
  FilterWarmStart::FilterWarmStart(const FilterWarmStart & other)
    :
    CoinWarmStartPrimalDual(other),
    CoinWarmStartBasis(other),
    lwsSize_(other.lwsSize_),
    lwsArray_(NULL),
    empty_(other.empty_)
  {
    assert(lwsSize_ > 0 || !lwsArray_);
    if (lwsSize_ > 0) {
      lwsArray_ = new fint[lwsSize_];
      assert(other.lwsArray_);
      CoinCopyN(other.lwsArray_, lwsSize_, lwsArray_);
    }
    for (int i = 0 ; i < 14 ; i ++) {
      istat_[i] = other.istat_[i];
    }
  }

  FilterWarmStart::~FilterWarmStart()
  {
    delete [] lwsArray_;
  }

  CoinWarmStartDiff *
  FilterWarmStart::generateDiff(const CoinWarmStart * const oldOne) const
  {
    const FilterWarmStart * old =
      dynamic_cast<const FilterWarmStart*> (oldOne);
    assert(old);

    CoinWarmStartDiff * diff = CoinWarmStartPrimalDual::generateDiff(old);

    CoinWarmStartPrimalDualDiff * pdDiff =
      dynamic_cast<CoinWarmStartPrimalDualDiff*>(diff);

    FilterWarmStartDiff* retval =
      new FilterWarmStartDiff(pdDiff, lwsSize_);
    delete diff;

    for (fint i = 0 ; i < lwsSize_ ; i++) {
      if (lwsArray_[i] != old->lwsArray_[i]) {
	retval->differences.push_back(FilterWarmStartDiff::OneDiff(i, lwsArray_[i] - old->lwsArray_[i]));
      }
    }

    retval->differences.resize(retval->differences.size());

    for (int i = 0 ; i < 14 ; i++) {
      retval->istat_[i] = istat_[i];
    }

    return retval;
  }

  void
  FilterWarmStart::applyDiff(const CoinWarmStartDiff * diff)
  {
    const FilterWarmStartDiff * diffF =
      dynamic_cast<const FilterWarmStartDiff *>(diff);
    assert(diffF);
    CoinWarmStartPrimalDual::applyDiff(diffF);

    fint end = static_cast<fint>(diffF->differences.size());
    for (fint i = 0 ; i < end ; i++) {
      lwsArray_[diffF->differences[i].first] += diffF->differences[i].second;
    }

    for (int i = 0 ; i < 14 ; i++)
      istat_[i] = diffF->istat_[i];
  }

  void
  FilterWarmStart::flushPoint()
  {
    delete [] lwsArray_;
  }

  FilterWarmStartDiff::FilterWarmStartDiff(CoinWarmStartPrimalDualDiff * diff,
					   fint capacity)
      :
      CoinWarmStartPrimalDualDiff()
  {
    CoinWarmStartPrimalDualDiff::swap(*diff);
    differences.reserve(capacity);
  }

  void
  FilterWarmStartDiff::flushPoint()
  {
    CoinWarmStartPrimalDualDiff::clear();
    differences.clear();
  }

  FilterTypes::fint
  FilterWarmStart::def_istat[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


} /* End namespace Bonmin */

