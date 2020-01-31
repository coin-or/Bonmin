// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/23/2007

#include "BonBabInfos.hpp"

namespace Bonmin
{
  /** Default constructor.*/
  BabInfo::BabInfo(int type):
      AuxInfo(type),
      babPtr_(NULL)
  {}

  /** Constructor from OsiBabSolver.*/
  BabInfo::BabInfo(const OsiBabSolver &other):
      AuxInfo(other),
      babPtr_(NULL)
  {}

  /** Copy constructor.*/
  BabInfo::BabInfo(const BabInfo &other):
      AuxInfo(other),
      babPtr_(other.babPtr_)
  {}

  /** Destructor.*/
  BabInfo::~BabInfo()
  {}

  /** Virtual copy constructor.*/
  OsiAuxInfo *
  BabInfo::clone() const
  {
    return new BabInfo(*this);
  }
}/* end namespace Bonmin*/

