// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
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
    CoinWarmStartBasis(),
    xSize_(xSize),
    xArray_(NULL),
    lamSize_(lamSize),
    lamArray_(NULL),
    lwsSize_(lwsSize),
    lwsArray_(NULL)
  {
    DBG_ASSERT(xSize > 0 || !xArray);
    if (xSize_ > 0){
      xArray_ = new real[xSize];
      DBG_ASSERT(xArray);
      CoinCopyN(xArray, xSize, xArray_);
    }
    DBG_ASSERT(lamSize > 0 || !lamArray);
    if (lamSize_ > 0){
      lamArray_ = new real[lamSize];
      DBG_ASSERT(lamArray);
      CoinCopyN(lamArray, lamSize, lamArray_);
    }
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
    CoinWarmStartBasis(other),
    xSize_(other.xSize_),
    xArray_(NULL),
    lamSize_(other.lamSize_),
    lamArray_(NULL),
    lwsSize_(other.lwsSize_),
    lwsArray_(NULL)
  {
    DBG_ASSERT(other.xSize_ > 0 || !other.xArray_);
    if (xSize_ > 0){
      xArray_ = new real[xSize_];
      DBG_ASSERT(other.xArray_);
      CoinCopyN(other.xArray_, xSize_, xArray_);
    }
    DBG_ASSERT(lamSize_ > 0 || !lamArray_);
    if (lamSize_ > 0){
      lamArray_ = new real[lamSize_];
      DBG_ASSERT(other.lamArray_);
      CoinCopyN(other.lamArray_, lamSize_, lamArray_);
    }
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
    delete [] xArray_;
    delete [] lamArray_;
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
  delete [] xArray_;
  delete [] lamArray_;
  delete [] lwsArray_;

  xArray_ = NULL;
  lamArray_ = NULL;
  lwsArray_ = NULL;
}

FilterSolver::fint 
BqpdWarmStart::def_istat[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


} /* End namespace Bonmin */
