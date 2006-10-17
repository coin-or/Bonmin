// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 02/15/2006


#include "BonIpoptWarmStart.hpp"
#include "CoinHelperFunctions.hpp"



#include "BonTMINLP2TNLP.hpp"
#include "BonIpoptInteriorWarmStarter.hpp"


namespace Bonmin {
/// Default constructor
IpoptWarmStart::IpoptWarmStart
(bool empty, int numvars, int numcont):
    CoinWarmStartBasis(),
    values_(),
    tempValues_(NULL),
    warm_starter_(NULL),
    empty_(empty)
{
  setSize(numvars,numcont);
}
/// Usefull constructor
IpoptWarmStart::IpoptWarmStart(const Ipopt::SmartPtr<TMINLP2TNLP> tnlp,
    SmartPtr<IpoptInteriorWarmStarter> warm_starter):
    values_(),
    tempValues_(NULL),
    warm_starter_(warm_starter),
    empty_(false)
{

  int numcols = tnlp->num_variables();
  int numrows = tnlp->num_constraints();
  setSize(numcols,numrows);
  values_.reserve(numcols+numrows);
  // For now, we keep all this in here, but we probably want to remove
  // it when we are happy with the new warmstarter (AW)
  double epsilon = 1e-05;//ipopt.getPushFact();
  const double * primals = tnlp->x_sol();
  const double * duals = tnlp->duals_sol();
  const double * colLo = tnlp->x_l();
  const double * colUp = tnlp->x_u();
  for(int i = 0 ; i < numcols ; i++) {
    if(primals[i] - colLo[i] < epsilon) {
      setStructStatus(i, atLowerBound);
      if(fabs(duals[i + numrows]) > epsilon) {
        values_.insert(i + numcols + numrows ,duals[i + numrows]);
      }

//       assert(duals[i +numrows + numcols] <= 100 * epsilon
// 	     ||
// 	     colUp[i] - colLo[i] <= epsilon
// 	       ||
// 	     colUp[i] - colLo[i] > 1e50);

    }
    else if( colUp[i] - primals[i] < epsilon) {
      setStructStatus(i, atUpperBound);
      if(fabs(duals[i + numrows + numcols]) > epsilon) {
        values_.insert(i + 2 * numcols + numrows ,duals[i + numrows + numcols]);
      }
//       assert(duals[i + numrows] <= epsilon
// 	     ||
// 	     colUp[i] - colLo[i] <= epsilon
// 	     ||
// 	     colUp[i] - colLo[i] > 1e100);
    }
    else {
      setStructStatus(i, basic);
//       assert((duals[i + numrows] <= epsilon && duals[i+ numrows +numcols] <= epsilon)
// 	     ||
// 	     colUp[i] - colLo[i] <= epsilon
// 	     ||
// 	     colUp[i] - colLo[i] > 1e100);

      values_.insert(i , primals[i]);
    }
  }

  // int i2 = 2*numcols;
  for(int i = 0 ; i < numrows ; i++) {
    if(fabs(duals[i])> epsilon) {
      values_.insert(i + numcols,duals[i]);
      setArtifStatus(i, basic);
    }
    else {
      setArtifStatus(i, atLowerBound);
    }

  }
  values_.sortIncrIndex();
}

/// Copy constructor
IpoptWarmStart::IpoptWarmStart( const IpoptWarmStart &other, bool ownValues):
    CoinWarmStartBasis(other),
    values_(other.values_),
    tempValues_(other.tempValues_),
    warm_starter_(NULL),//(other.warm_starter_),
    empty_(other.empty_)
{
  //  if(ownValues_ && other.values_ != NULL)
}


CoinWarmStartDiff*
IpoptWarmStart::generateDiff(const CoinWarmStart *const oldCWS) const
{
  CoinWarmStartDiff * diff = CoinWarmStartBasis::generateDiff(oldCWS);
  CoinWarmStartBasisDiff * basisDiff =
    dynamic_cast<CoinWarmStartBasisDiff *>(diff);


  CoinWarmStartDiff* retval =
    new IpoptWarmStartDiff(basisDiff, values_, NULL);//warm_starter_);
  delete diff;
  return retval;
}


void
IpoptWarmStart::applyDiff (const CoinWarmStartDiff *const cwsdDiff)
{
  CoinWarmStartBasis::applyDiff(cwsdDiff);
  IpoptWarmStartDiff const * const ipoptDiff =
    dynamic_cast<IpoptWarmStartDiff const * const > (cwsdDiff);

  tempValues_ = ipoptDiff->diffValues_;
  //  ownValues_ = 0;
  warm_starter_ = ipoptDiff->warm_starter();
}

IpoptWarmStart::~IpoptWarmStart()
{}

void
IpoptWarmStart::flushPoint()
{
  if(values_.getNumElements() > 0)
    values_.clear();
}
void

IpoptWarmStartDiff::flushPoint()
{
  if(diffValues_) {
    delete diffValues_;
    diffValues_ = NULL;
  }
}
}
