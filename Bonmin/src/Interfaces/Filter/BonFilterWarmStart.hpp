// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 11/21/2006




#ifndef BonFilterWarmStart_H
#define BonFilterWarmStart_H

#include "CoinWarmStartBasis.hpp"

namespace Bonmin{

  class BonFilterWarmStart : public CoinWarmStartBasis
  {
  };


  class BonFilterWarmStartDiff : public CoinWarmStartBasisDiff
  {
  };

} /* end namespace Bonmin */
#endif

