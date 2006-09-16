#include "CoinHelperFunctions.hpp"
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

    if (lower_bound_ >= upper_bound() + cutoffIncr) {
	return BCP_DoNotBranch_Fathomed;
    }

    if (numNlpFailed_ >= par.entry(BM_par::NumNlpFailureMax)) {
	printf("WARNING! Too many (%i) NLP failures in a row. Abandoning node.",
	       numNlpFailed_);
	return BCP_DoNotBranch_Fathomed;
    }

    if (par.entry(BM_par::BranchOnSos) && sos.size() > 0) {
	double bestval = 2.0;
	int bestsos = -1;
	int bestsplit = -1;
	int bestprio = 0;
	for (size_t ind = 0; ind < sos.size(); ++ind) {
	    const BmSosInfo& sosi = *sos[ind];
	    const int prio = sosi.priority;
	    const int type = sosi.type;
	    if ((type == 0 && par.entry(BM_par::BranchOnSos) == 2) ||
		(type == 1 && par.entry(BM_par::BranchOnSos) == 1))
		continue;
	    if (bestsos >= 0 && prio != bestprio)
		break;
	    const int first = sosi.first;
	    const int last = sosi.last;
	    const int *indices = sosi.indices;
	    // First count how many nonzeros are there
	    int cnt = 0;
	    for (int j = first; j < last; ++j) {
		const double prim_i = primal_solution_[indices[j]];
		if (prim_i > intTol && prim_i < 1-intTol) {
		    ++cnt;
		}
	    }
	    // For type1 there can be at most 1, for type2 at most 2 nonzeros
	    // Find the next sos if this is satisfied.
	    if (cnt <= 1+type)
		continue;

	    // Find where does the sos goes over 0.5
	    int half = -10;
	    double val = 0.0;
	    double primal = 0.0;
	    for (half = first; half < last; ++half) {
		val = primal_solution_[indices[half]];
		primal += val;
		if (primal > 0.5) {
		    break;
		}
	    }
	    /* Now [first,half) is < 0.5 and [first,half] is > 0.5.
	       If the latter is closer to 0.5 then bump up half by one. */
	    if (primal - 0.5 < 0.5 - (primal - val)) {
		++half;
		val = primal - 0.5;
	    } else {
		val = 0.5 - (primal - val);
	    }

	    /* Now [first,half) gets as close to 0.5 as possible, and val
	       contains its distance from 0.5 */

	    /* pick this over the prev best choice only if val is smaller than
	       the it was for the prev best. */
	    if (val >= bestval)
		continue;

	    /* The branches will be as follows:
	       type1 : [first,half) and [half,last)
	       type2 : [first,half) and (half,last)
	       So if half==first then it must be incremented,
	       and also for type 2 if half==last-1 then it must be decremented.
	    */
	    if (half == first)
		++half;
	    if (type == 1 && half == last-1)
		--half;

	    bestval = val;
	    bestsos = ind;
	    bestsplit = half;
	    bestprio = prio;
	}

	if (bestsos >= 0) {
	    // OK, we can branch on an SOS
	    // This vector will contain one entry only, but we must pass in a
	    // vector
	    BCP_vec<BCP_var*> new_vars;
	    // This vector contains which vars have their ounds forcibly
	    // changed. It's going to be the newly added branching var, hence
	    // it's position is -1
	    BCP_vec<int> fvp(1, -1);
	    // This vector contains the new lb/ub pairs for all children. So
	    // it's going to be {0.0, 0.0, 1.0, 1.0}
	    BCP_vec<double> fvb(4, 0.0);
	    fvb[2] = fvb[3] = 1.0;
	    new_vars.push_back(new BM_branching_var(bestsos, bestsplit));
	    cands.push_back(new BCP_lp_branching_object(2, // num of children
							&new_vars,
							NULL, // no new cuts
							&fvp,NULL,&fvb,NULL,
							NULL,NULL,NULL,NULL));
	    if (par.entry(BM_par::PrintBranchingInfo)) {
		printf("\
BM_lp: At node %i branching on SOS%i set %i  split at: %i with val: %f\n",
		       current_index(), sos[bestsos]->type+1,
		       bestsos, bestsplit, bestval);
	    }
	    return BCP_DoBranch;;
	}
    }

    /* So we couldn't do an sos branching. If it's not pure B&B then we might
       as well do the built-in branching */
    if (! par.entry(BM_par::PureBranchAndBound)) {
	return BCP_lp_user::select_branching_candidates(lpres, vars, cuts,
							local_var_pool,
							local_cut_pool, cands);
    }

    /* If PURE_BB then we have the NLP optimum in "primal_solution_". 
       Try to find a good branching candidate there. */
    const int numVars = vars.size();
    double* frac = new double[numVars];
    for (int i = 0; i < numVars; ++i) {
	if (vars[i]->var_type() == BCP_ContinuousVar) {
	    frac[i] = 0.0;
	    continue;
	}
	const double psol = CoinMin(CoinMax(vars[i]->lb(),primal_solution_[i]),
				    vars[i]->ub());
	const double f = psol - floor(psol);
	frac[i] = CoinMin(f, 1.0-f);
    }

    // if (par.entry(BM_par::StrongBranchCandidateMinFrac) > 0) {}
    
    int besti = -1;
    double bestval = -1e200;
    for (int i = 0; i < numVars; ++i) {
	if (vars[i]->var_type() == BCP_ContinuousVar)
	    continue;
	const double psol = CoinMin(CoinMax(vars[i]->lb(),primal_solution_[i]),
				    vars[i]->ub());
	const double frac = psol - floor(psol);
	const double dist = CoinMin(frac, 1.0-frac);
	
	if (dist < intTol)
	    continue;
	if (par.entry(BM_par::CombinedDistanceAndPriority)) {
	    if (dist * branching_priority_[i] > bestval) {
		bestval = dist * branching_priority_[i];
		besti = i;
	    }
	}
	else {
	    if (branching_priority_[i] > bestval) {
		bestval = branching_priority_[i];
		besti = i;
	    }
	}
    }
    if (besti == -1) {
	return BCP_DoNotBranch_Fathomed;
    }
    double psol = CoinMin(CoinMax(vars[besti]->lb(), primal_solution_[besti]),
			  vars[besti]->ub());
    if (par.entry(BM_par::PrintBranchingInfo)) {
	printf("BM_lp: branching on variable %i   value: %f\n", besti, psol);
    }
    BCP_vec<int> select_pos(1, besti);
    append_branching_vars(primal_solution_, vars, select_pos, cands);
    return BCP_DoBranch;
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
