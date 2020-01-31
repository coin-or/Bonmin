// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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

using namespace Ipopt;

namespace Bonmin
{
/// Default constructor
  IpoptWarmStart::IpoptWarmStart
  (bool empty, int numvars, int numcont):
      CoinWarmStartPrimalDual(),
      CoinWarmStartBasis(),
      warm_starter_(NULL),
      empty_(empty)
  {
    setSize(numvars,numcont);
  }
/// Usefull constructor
  IpoptWarmStart::IpoptWarmStart(const Ipopt::SmartPtr<TMINLP2TNLP> tnlp,
      SmartPtr<IpoptInteriorWarmStarter> warm_starter):
      CoinWarmStartPrimalDual(tnlp->num_variables(),
			      2*tnlp->num_variables()+tnlp->num_constraints(),
			      tnlp->x_sol(), tnlp->duals_sol() ),
      CoinWarmStartBasis(),
      warm_starter_(warm_starter),
      empty_(false)
  {
    int numcols = tnlp->num_variables();
    int numrows = tnlp->num_constraints();
    // initialize the size in the CoinWarmStartBasis part
    setSize(numcols,numrows);
  }

/// Another usefull constructor, stores the passed point
IpoptWarmStart::IpoptWarmStart(int primal_size, int dual_size,
                   const double * primal, const double * dual):
          CoinWarmStartPrimalDual(primal_size, dual_size, primal, dual),
          CoinWarmStartBasis(),
          warm_starter_(NULL), 
          empty_(false)
{
   setSize(primal_size, dual_size - 2* primal_size);
}
/// Copy constructor
  IpoptWarmStart::IpoptWarmStart( const IpoptWarmStart &other, bool ownValues):
    CoinWarmStartPrimalDual(other),
    CoinWarmStartBasis(other),
    warm_starter_(NULL /*other.warm_starter_*/),
    empty_(other.empty_)
  {
    //  if(ownValues_ && other.values_ != NULL)
  }

/// A constructor from a CoinWarmStartPrimalDual
  IpoptWarmStart::IpoptWarmStart(const CoinWarmStartPrimalDual& pdws) :
    CoinWarmStartPrimalDual(pdws),
    CoinWarmStartBasis(),
    warm_starter_(NULL),
    empty_(false)
  {   
  }
  
  CoinWarmStartDiff*
  IpoptWarmStart::generateDiff(const CoinWarmStart *const oldCWS) const
  {
    const IpoptWarmStart * const ws =
      dynamic_cast< const IpoptWarmStart * const > (oldCWS);
    DBG_ASSERT(ws);

    CoinWarmStartDiff * diff = CoinWarmStartPrimalDual::generateDiff(ws);

    CoinWarmStartPrimalDualDiff * pdDiff =
      dynamic_cast<CoinWarmStartPrimalDualDiff*>(diff);

    CoinWarmStartDiff* retval =
      new IpoptWarmStartDiff(pdDiff, NULL);//warm_starter_);
    delete diff;

    return retval;
  }


  void
  IpoptWarmStart::applyDiff (const CoinWarmStartDiff *const cwsdDiff)
  {
    IpoptWarmStartDiff const * const ipoptDiff =
      dynamic_cast<IpoptWarmStartDiff const * const > (cwsdDiff);
    DBG_ASSERT(ipoptDiff);
    CoinWarmStartPrimalDual::applyDiff(ipoptDiff);
    warm_starter_ = ipoptDiff->warm_starter();
  }

  IpoptWarmStart::~IpoptWarmStart()
  {}

  void
  IpoptWarmStart::flushPoint()
  {
    CoinWarmStartPrimalDual::clear();
  }

  void
  IpoptWarmStartDiff::flushPoint()
  {
    CoinWarmStartPrimalDualDiff::clear();
  }
}
