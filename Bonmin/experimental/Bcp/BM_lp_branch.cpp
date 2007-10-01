// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "CoinHelperFunctions.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_lp.hpp"
#include "BM.hpp"
#include "BonChooseVariable.hpp"
#include "BonCurvBranchingSolver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"
#include "BonOsiTMINLPInterface.hpp"

//#############################################################################

BCP_branching_decision
BM_lp::select_branching_candidates(const BCP_lp_result& lpres,
                                   const BCP_vec<BCP_var*>& vars,
                                   const BCP_vec<BCP_cut*>& cuts,
                                   const BCP_lp_var_pool& local_var_pool,
                                   const BCP_lp_cut_pool& local_cut_pool,
                                   BCP_vec<BCP_lp_branching_object*>& cands,
				   bool force_branch)
{
  Bonmin::OsiTMINLPInterface* nlp =
    dynamic_cast<Bonmin::OsiTMINLPInterface*>(getLpProblemPointer()->lp_solver);
  // If we are doing pure B&B then we have an nlp, and then we check for
  // consecutive failures.
  if (nlp) {
    if (lpres.termcode() & BCP_Abandoned) {
      if (nlp->isIterationLimitReached()) {
	printf("\
BM_lp: At node %i : WARNING: nlp reached iter limit. Will force branching\n",
	       current_index());
      } else {
	printf("\
BM_lp: At node %i : WARNING: nlp is abandoned. Will force branching\n",
	       current_index());
      }
      // nlp failed
      nlp->forceBranchable();
      numNlpFailed_++;
      if (numNlpFailed_ >= par.entry(BM_par::NumNlpFailureMax)) {
	printf("WARNING! Too many (%i) NLP failures in a row. Abandoning node.",
	       numNlpFailed_);
	return BCP_DoNotBranch_Fathomed;
      }
    } else {
      numNlpFailed_ = 0;
    }
  }

  OsiBranchingInformation brInfo(nlp, false, true);
  brInfo.cutoff_ = upper_bound() + get_param(BCP_lp_par::Granularity);
  brInfo.integerTolerance_ = integerTolerance_;
  brInfo.timeRemaining_ = get_param(BCP_lp_par::MaxRunTime) - CoinCpuTime();
  brInfo.numberSolutions_ = 0; /*FIXME*/
  brInfo.numberBranchingSolutions_ = 0; /*FIXME numBranchingSolutions_;*/
  brInfo.depth_ = current_level();

  BCP_branching_decision brDecision;
  if (bonmin_.getAlgorithm() == 0) {
    /* if pure B&B */
    brDecision = bbBranch(brInfo, cands);
  } else {
    brDecision = hybridBranch();
  }

  return brDecision;
}

//-----------------------------------------------------------------------------

BCP_branching_decision
BM_lp::hybridBranch()
{
    // FIXME: most of the pureBB stuff should work here.
    throw BCP_fatal_error("BM_lp: FIXME: make hybrid work...");
}

//-----------------------------------------------------------------------------

BCP_branching_decision
BM_lp::bbBranch(OsiBranchingInformation& brInfo,
		BCP_vec<BCP_lp_branching_object*>& cands)
{
  OsiSolverInterface* osi = getLpProblemPointer()->lp_solver;
  Bonmin::OsiTMINLPInterface* nlp =
    dynamic_cast<Bonmin::OsiTMINLPInterface*>(osi);
  assert(nlp);

  nlp->getDblParam(OsiPrimalTolerance, brInfo.integerTolerance_);
    
  BCP_branching_decision retCode;
  OsiBranchingObject* brObj = NULL;

//   static int cnt = 0;
//   printf("cnt = %i\n", cnt);
//   ++cnt;

  const int numCols = nlp->getNumCols();
  double* clb_old = new double[numCols];
  double* cub_old = new double[numCols];
  CoinDisjointCopyN(nlp->getColLower(), numCols, clb_old);
  CoinDisjointCopyN(nlp->getColUpper(), numCols, cub_old);

  Ipopt::SmartPtr<Ipopt::OptionsList> options = bonmin_.options();
  int numSB = 0;
  const bool sbIsSet =
    options->GetIntegerValue("number_strong_branch",numSB,"bonmin.");
  int numSBroot = 0;
  const bool sbRootIsSet =
    options->GetIntegerValue("number_strong_branch_root",numSBroot,"bonmin.");

  if (sbRootIsSet && brInfo.depth_ == 0) {
    bonmin_.branchingMethod()->setNumberStrong(numSBroot);
  } else {
    bonmin_.branchingMethod()->setNumberStrong(numSB);
  }

  const int brResult =
    try_to_branch(brInfo, nlp, bonmin_.branchingMethod(), brObj, true);

#if 0
  if (choose->goodSolution()) {
    /* yipee! a feasible solution! Check that it's really */
    const double* sol = choose->goodSolution();
    BM_solution* bmsol = new BM_solution;
    //Just copy the solution stored in solver to sol
    double ptol = integerTolerance_;
    for (int i = 0 ; i < numCols ; i++) {
      if (fabs(sol[i]) > ptol)
	bmsol->add_entry(i, sol[i]); 
    }
    bmsol->setObjective(choose->goodObjectiveValue());
    choose->clearGoodSolution();
    send_feasible_solution(bmsol);
    delete bmsol;
  }
#endif

  switch (brResult) {
  case -2:
    // when doing strong branching a candidate has proved that the
    // problem is infeasible
    retCode = BCP_DoNotBranch_Fathomed;
    break;
  case -1:
    // OsiChooseVariable::chooseVariable() returned 2, 3, or 4
    if (!brObj) {
      // just go back and resolve
      retCode = BCP_DoNotBranch;
    } else {
      // otherwise might as well branch. The forced variable is
      // unlikely to jump up one more (though who knows...)
      retCode = BCP_DoBranch;
    }
    break;
  case 0:
    if (!brObj) {
      // nothing got fixed, yet couldn't find something to branch on
      throw BCP_fatal_error("BM: Couldn't branch!\n");
    }
    // we've got a branching object
    retCode = BCP_DoBranch;
    break;
  default:
    throw BCP_fatal_error("\
BM: BCP_lp_user::try_to_branch returned with unknown return code.\n");
  }

  bool distributedStrongBranching =
    (retCode == BCP_DoBranch) && (par.entry(BM_par::FullStrongBranch) > 0);

  if (brResult < 0) {
    // If there were some fixings (brResult < 0) then propagate them
    // where needed

    // FIXME: This is not nice. Meddling w/ BCP internal data. The BCP
    // user interface should provide a way to change bounds regardless
    // whether branching is asked for or not.
    const double* clb = nlp->getColLower();
    const double* cub = nlp->getColUpper();

    BCP_vec<BCP_var*>& vars = getLpProblemPointer()->node->vars;
    if (distributedStrongBranching) {
      retCode = retCode;
    }
    for (int i = 0; i < numCols; ++i) {
      if (clb_old[i] != clb[i] || cub_old[i] != cub[i]) {
	vars[i]->change_bounds(clb[i], cub[i]);
	nlp->setColBounds(i, clb[i], cub[i]);
      }
    }
  }

  if (retCode == BCP_DoBranch) {
    // all possibilities are 2-way branches
    int order[2] = {0, 1};
    // Now interpret the result (at this point we must have a brObj
    OsiIntegerBranchingObject* intBrObj =
      dynamic_cast<OsiIntegerBranchingObject*>(brObj);
    if (intBrObj) {
      if (intBrObj->firstBranch() == 1) {
	order[0] = 1;
	order[1] = 0;
      }
      BCP_lp_integer_branching_object o(intBrObj);
      cands.push_back(new BCP_lp_branching_object(o, order));
      if (par.entry(BM_par::PrintBranchingInfo)) {
	printf("BM_lp: branching on variable %i   value: %f\n",
	       intBrObj->originalObject()->columnNumber(),
	       intBrObj->value());
      }
    }
    OsiSOSBranchingObject* sosBrObj =
      dynamic_cast<OsiSOSBranchingObject*>(brObj);
    if (sosBrObj) {
      if (sosBrObj->firstBranch() == 1) {
	order[0] = 1;
	order[1] = 0;
      }
      BCP_lp_sos_branching_object o(sosBrObj);
      cands.push_back(new BCP_lp_branching_object(nlp, o, order));
      if (par.entry(BM_par::PrintBranchingInfo)) {
	printf("BM_lp: branching on SOS %i   value: %f\n",
	       sosBrObj->originalObject()->columnNumber(),
	       sosBrObj->value());
      }
    }
  }

  delete brObj;
  delete[] clb_old;
  delete[] cub_old;
  return retCode;
}

/****************************************************************************/

void
BM_lp::set_user_data_for_children(BCP_presolved_lp_brobj* best, 
                                  const int selected)
{
    BM_node* data = NULL;
    data = new BM_node;
    data->numNlpFailed_ = numNlpFailed_;
    best->user_data()[0] = data;
    data = new BM_node;
    data->numNlpFailed_ = numNlpFailed_;
    best->user_data()[1] = data;
}

//#############################################################################

void
BM_lp::process_message(BCP_buffer& buf)
{
}

/*
 Ipopt::TMINLP* model = nlp.model();
 model->get_nlp_info to get nnz_h_lag (pp 16)
alloc that much space for irow, jcol and values  (passed into eval_h, pp22)
call eval_h *twice*
first: values=NULL to get the sparsity structure
second: irow and jcol are NULL and values will be filled
the sparsity format is in MA27 format (pp 39 -- but it's symmetric here!)
(multiple triplets can point to the same location! -- add them up!)

the diagonal of that matrix is the 1-dim quadratic term for each var when the
others are fixed. Better be nonneg for convex problems)
*/
