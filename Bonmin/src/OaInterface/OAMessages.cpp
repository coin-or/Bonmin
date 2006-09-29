// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005


#include "OAMessages.hpp"

namespace Bonmin{

OaMessages::OaMessages():
    CoinMessages(DUMMY_END)
{
  strcpy(source_,"OA");
  addMessage(FEASIBLE_NLP,CoinOneMessage( 1, 2,"Solved NLP in %d iterations, found a feasible solution of value %f."));
  addMessage(INFEASIBLE_NLP, CoinOneMessage(2,2,"Solved NLP in %d iterations, problem is infeasible in subspace."));
  addMessage(UPDATE_UB, CoinOneMessage(3,1,"New best feasible of %g found after %g sec."));
  addMessage(SOLVED_LOCAL_SEARCH, CoinOneMessage(4,2,"Local search solved to optimality in %d nodes and %d lp iterations."));
  addMessage(LOCAL_SEARCH_ABORT, CoinOneMessage(5,2,"Local search aborted : %d nodes and %d lp iterations."));
  addMessage(UPDATE_LB, CoinOneMessage(6,2,"Updating lower bound to %g elapsed time %g sec"));
  addMessage(ABORT,CoinOneMessage(7,1,"Oa aborted on %s limit, time spent %g"));
  addMessage(OASUCCESS, CoinOneMessage(8,1,"Oa converged in %g seconds"));
  addMessage(LP_ERROR,CoinOneMessage(9,2,"Error of LP approximation %g"));
  addMessage(PERIODIC_MSG, CoinOneMessage(10,1,"After %7.1f seconds, upper bound %10g, lower bound %10g"));
}

}//end namespace Bonmin
