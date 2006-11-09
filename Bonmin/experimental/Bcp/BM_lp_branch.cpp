#include "CoinHelperFunctions.hpp"
#include "BCP_lp_node.hpp"
#include "BM.hpp"

//#############################################################################

BCP_branching_decision
BM_lp::select_branching_candidates(const BCP_lp_result& lpres,
                                   const BCP_vec<BCP_var*>& vars,
                                   const BCP_vec<BCP_cut*>& cuts,
                                   const BCP_lp_var_pool& local_var_pool,
                                   const BCP_lp_cut_pool& local_cut_pool,
                                   BCP_vec<BCP_lp_branching_object*>& cands)
{
    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
    double intTol;
    double cutoffIncr;
    options->GetNumericValue("integer_tolerance", intTol, "bonmin.");
    options->GetNumericValue("cutoff_decr", cutoffIncr, "bonmin.");

    const double objLimit = upper_bound() + cutoffIncr;

    if (lower_bound_ >= objLimit) {
	return BCP_DoNotBranch_Fathomed;
    }

    if (numNlpFailed_ >= par.entry(BM_par::NumNlpFailureMax)) {
	printf("WARNING! Too many (%i) NLP failures in a row. Abandoning node.",
	       numNlpFailed_);
	return BCP_DoNotBranch_Fathomed;
    }

    // FIXME: for now let's work on pure B&B, though most of this will work
    // FIXME: for general hybrid algo
    if (par.entry(BM_par::PureBranchAndBound)) {
	OsiBranchingObject* brObj = NULL;
	OsiBranchingInformation brInfo(&nlp, false);

	brInfo.cutoff_ = objLimit;
	brInfo.integerTolerance_ = intTol;
	brInfo.timeRemaining_ =
	    get_param(BCP_lp_par::MaxRunTime) - CoinCpuTime();

	brInfo.numberSolutions_ = 0; /*FIXME*/
	brInfo.numberBranchingSolutions_ = 0; /*FIXME numBranchingSolutions_;*/
	brInfo.depth_ = current_level();

	const int numCols = nlp.getNumCols();
	double* clb_old = new double[numCols];
	double* cub_old = new double[numCols];
	CoinDisjointCopyN(nlp.getColLower(), numCols, clb_old);
	CoinDisjointCopyN(nlp.getColUpper(), numCols, cub_old);
	
	OsiChooseStrong choose(&nlp);
	choose.setNumberBeforeTrusted(5); // the default in Cbc
	choose.setNumberStrong(5); // the default in Cbc
	/** Pseudo Shadow Price mode
	    0 - off
	    1 - use and multiply by strong info
	    2 - use 
	*/
	choose.setShadowPriceMode(0);

	const int brResult = try_to_branch(brInfo, &nlp, &choose, brObj, true);
#if 0
	/* FIXME:before doing anything check if we have found a new solution */
	if (choose->goodSolution()
	    &&model->problemFeasibility()->feasible(model,-1)>=0) {
	    // yes
	    double objValue = choose->goodObjectiveValue();
	    model->setBestSolution(CBC_STRONGSOL,
				   objValue,
				   choose->goodSolution()) ;
	    model->setLastHeuristic(NULL);
	    model->incrementUsed(choose->goodSolution());
	    choose->clearGoodSolution();
	}
#endif
	switch (brResult) {
	case -2:
	    // when doing strong branching a candidate has proved that the
	    // problem is infeasible
	    delete[] clb_old;
	    delete[] cub_old;
	    return BCP_DoNotBranch_Fathomed;
	case -1:
	    // OsiChooseVariable::chooseVariable() returned 2, 3, or 4
	    if (!brObj) {
		// just go back and resolve
		delete[] clb_old;
		delete[] cub_old;
		return BCP_DoNotBranch;
	    }
	    // otherwise might as well branch. The forced variable is
	    // unlikely to jump up one more (though who knows...)
	    break;
	case 0:
	    if (!brObj) {
		// nothing got fixed, yet couldn't find something to branch on
		throw BCP_fatal_error("BM: Couldn't branch!\n");
	    }
	    // we've got a branching object
	    break;
	default:
	    throw BCP_fatal_error("\
BM: BCP_lp_user::try_to_branch returned with unknown return code.\n");
	}

	// If there were some fixings (brResult < 0) then propagate them where
	// needed
	if (brResult < 0) {
	    const double* clb = nlp.getColLower();
	    const double* cub = nlp.getColUpper();
	    BCP_lp_prob* p = getLpProblemPointer();
	    BCP_vec<BCP_var*>& vars = p->node->vars;
	    OsiSolverInterface* lp = p->lp_solver;
	    for (int i = 0; i < numCols; ++i) {
		if (clb_old[i] != clb[i] || cub_old[i] != cub[i]) {
		    vars[i]->change_bounds(clb[i], cub[i]);
		    lp->setColBounds(i, clb[i], cub[i]);
		}
	    }
	}
	delete[] clb_old;
	delete[] cub_old;

	// Now interpret the result (at this point we must have a brObj
	OsiIntegerBranchingObject* intBrObj =
	    dynamic_cast<OsiIntegerBranchingObject*>(brObj);
	if (intBrObj) {
	    BCP_lp_integer_branching_object o(intBrObj);
	    cands.push_back(new BCP_lp_branching_object(o));
	    if (par.entry(BM_par::PrintBranchingInfo)) {
		printf("BM_lp: branching on variable %i   value: %f\n",
		       intBrObj->originalObject()->columnNumber(),
		       intBrObj->value());
	    }
	}
	OsiSOSBranchingObject* sosBrObj =
	    dynamic_cast<OsiSOSBranchingObject*>(brObj);
	if (sosBrObj) {
	    BCP_lp_sos_branching_object o(sosBrObj);
	    cands.push_back(new BCP_lp_branching_object(&nlp, o));
	    if (par.entry(BM_par::PrintBranchingInfo)) {
		printf("BM_lp: branching on variable %i   value: %f\n",
		       sosBrObj->originalObject()->columnNumber(),
		       sosBrObj->value());
	    }
	}
	return BCP_DoBranch;

    } else { // not pure B&B
	throw BCP_fatal_error("BM_lp: FIXME: make hybrid work...");
    }

    return BCP_DoBranch; // to quiet gcc
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
