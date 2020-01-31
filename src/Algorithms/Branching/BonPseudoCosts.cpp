// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007

#include "BonPseudoCosts.hpp"


namespace Bonmin
{

  PseudoCosts::PseudoCosts():
    OsiPseudoCosts()
  {
  }

  /** Copy constructor.*/
  PseudoCosts::PseudoCosts(const PseudoCosts & rhs):
   OsiPseudoCosts(rhs)
  {
  }


  /** Assignment operator const version.*/
  PseudoCosts &
  PseudoCosts::operator=(const PseudoCosts&rhs)
  {
    if (this != &rhs) {
        OsiPseudoCosts::operator=(rhs);
    }
    return *this;
  }

}/* End Bonmin namespace.*/

