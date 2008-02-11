// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//
// Date : 02/10/2008

#include "BonCouenneInfo.hpp"

namespace Bonmin
{
  /** Default constructor.*/
  CouenneInfo::CouenneInfo(int type):
      BabInfo(type)
  {}

  /** Constructor from OsiBabSolver.*/
  CouenneInfo::CouenneInfo(const OsiBabSolver &other):
      BabInfo(other)
  {}

  /** Copy constructor.*/
  CouenneInfo::CouenneInfo(const CouenneInfo &other):
      BabInfo(other)
  {}

  /** Destructor.*/
  CouenneInfo::~CouenneInfo()
  {}

  /** Virtual copy constructor.*/
  OsiAuxInfo *
  CouenneInfo::clone() const
  {
    return new CouenneInfo(*this);
  }

  CouenneInfo::NlpSolution::NlpSolution(int n, const double* sol, double objval)
    :
    n_(n),
    objVal_(objval)
  {
    sol_ = new double[n];
    CoinCopyN(sol, n, sol_);
  }

  CouenneInfo::NlpSolution::~NlpSolution()
  {
    delete [] sol_;
  }
}/* end namespace Bonmin*/

