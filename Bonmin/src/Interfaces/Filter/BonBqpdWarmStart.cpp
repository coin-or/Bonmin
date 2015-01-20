// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                    based on BonFilterWarmStart.cpp
//
// Date : 2007-08-03

#include "BonBqpdWarmStart.hpp"


namespace Bonmin{

  BqpdWarmStart:: BqpdWarmStart(const fint xSize /*= 0*/,
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
    DBG_ASSERT(lwsSize > 0 || !lwsArray);
    if (lwsSize_ > 0){
      lwsArray_ = new fint[lwsSize];
      DBG_ASSERT(lwsArray);
      CoinCopyN(lwsArray, lwsSize, lwsArray_);
    }
    for(int i = 0 ; i < 14 ; i ++) {
      istat_[i] = istat[i];
    }
  }

  /* Copy constructor */
  BqpdWarmStart::BqpdWarmStart(const BqpdWarmStart & other)
    :
    CoinWarmStartPrimalDual(other),
    CoinWarmStartBasis(other),
    lwsSize_(other.lwsSize_),
    lwsArray_(NULL),
    empty_(other.empty_)
  {
    DBG_ASSERT(lwsSize_ > 0 || !lwsArray_);
    if (lwsSize_ > 0){
      lwsArray_ = new fint[lwsSize_];
      DBG_ASSERT(other.lwsArray_);
      CoinCopyN(other.lwsArray_, lwsSize_, lwsArray_);
    }
    for(int i = 0 ; i < 14 ; i ++) {
      istat_[i] = other.istat_[i];
    }
  }

  BqpdWarmStart::~BqpdWarmStart()
  {
    delete [] lwsArray_;
  }

  CoinWarmStartDiff *
  BqpdWarmStart::generateDiff(const CoinWarmStart * const oldOne) const
  {
    throw CoinError("Method not implemented",
		    "generateDiffs",
		    "BqpdWarmStart");
    return NULL;
  }


  void
  BqpdWarmStart::applyDiff(const CoinWarmStartDiff * diff){
    
    throw CoinError("Method not implemented",
		    "applyDiff",
		    "BqpdWarmStart");
  }

void
BqpdWarmStart::flushPoint()
{
  CoinWarmStartPrimalDual::clear();
  delete [] lwsArray_;

  lwsArray_ = NULL;
}

FilterSolver::fint 
BqpdWarmStart::def_istat[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


} /* End namespace Bonmin */
