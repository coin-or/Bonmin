// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005


#include "BonOAMessages.hpp"
#include <cstring>
#include "BonMsgUtils.hpp"

namespace Bonmin
{

  OaMessages::OaMessages():
      CoinMessages(DUMMY_END)
  {
    strcpy(source_,"OA");
    ADD_MSG(FEASIBLE_NLP, std_m, 2,"Solved NLP in %d iterations, found a feasible solution of value %f.");
    ADD_MSG(INFEASIBLE_NLP, std_m,2,"Solved NLP in %d iterations, problem is infeasible in subspace.");
    ADD_MSG(UPDATE_UB, std_m,1,"New best feasible of %g found after %g sec and %i iterations.");
    ADD_MSG(SOLVED_LOCAL_SEARCH, std_m,2,"Local search solved to optimality in %d nodes and %d lp iterations.");
    ADD_MSG(LOCAL_SEARCH_ABORT, std_m,2,"Local search aborted : %d nodes and %d lp iterations.");
    ADD_MSG(UPDATE_LB, std_m ,2,"Updating lower bound to %g elapsed time %g sec");
    ADD_MSG(ABORT,std_m,1,"%s aborted on %s limit, time spent %g has solution of value %g");
    ADD_MSG(OASUCCESS, std_m,1,"%s converged in %g seconds found solution of value %g (lower bound %g ).");
    ADD_MSG(OAABORT, std_m,1,"%s interupted after %g seconds found solution of value %g (lower bound %g ).");
    ADD_MSG(OA_STATS, std_m,1,"Performed %i iterations, explored %i branch-and-bound nodes in total");
    ADD_MSG(LP_ERROR,std_m,2,"Error of LP approximation %g");
    ADD_MSG(PERIODIC_MSG, std_m,1,"After %7.1f seconds, %i iterations upper bound %10g, lower bound %10g");
    ADD_MSG(FP_DISTANCE, std_m,2,"NLP solution is %g from MILP sol");
    ADD_MSG(FP_MILP_VAL, std_m,2,"MILP solution has value w.r.t original objective: %10g");
    ADD_MSG(FP_MAJOR_ITERATION, std_m,1,"Major iteration %i ub: %g");
    ADD_MSG(FP_MINOR_ITERATION, std_m,1,"Minor iteration %i ub: %g");
  }

}//end namespace Bonmin
