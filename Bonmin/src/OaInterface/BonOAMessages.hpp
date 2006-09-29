// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  07/15/2005

#ifndef OaMessages_H
#define OaMessages_H
#include "CoinMessage.hpp"

namespace Bonmin
{
  enum OA_Message{
    FEASIBLE_NLP,
    INFEASIBLE_NLP,
    UPDATE_UB,
    SOLVED_LOCAL_SEARCH,
    LOCAL_SEARCH_ABORT,
    UPDATE_LB,
    ABORT,
    OASUCCESS,
    LP_ERROR,
    PERIODIC_MSG,
    DUMMY_END
  };

  /** Output messages for Outer approximation cutting planes */
  class OaMessages : public CoinMessages
  {
  public:
    OaMessages();
  };

} //end namespace Bonmin
#endif
