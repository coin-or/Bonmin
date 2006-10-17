// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
//

#include "BM.hpp"

//#############################################################################

void
BM_lp::pack_feasible_solution(BCP_buffer& buf, const BCP_solution* sol)
{
    const BM_solution* bs = dynamic_cast<const BM_solution*>(sol);
    if (!bs) {
	throw BCP_fatal_error("Trying to pack non-BM_solution.\n");
    }
    buf.pack(bs->_objective);
    buf.pack(bs->_ind);
    buf.pack(bs->_values);
}

/****************************************************************************/

BCP_solution*
BM_tm::unpack_feasible_solution(BCP_buffer& buf)
{
    BM_solution* bs = new BM_solution;
    buf.unpack(bs->_objective);
    buf.unpack(bs->_ind);
    buf.unpack(bs->_values);
    return bs;
}

//#############################################################################

void
BM_tm::pack_module_data(BCP_buffer& buf, BCP_process_t ptype)
{
    // possible process types looked up in BCP_enum_process_t.hpp
    switch (ptype) {
    case BCP_ProcessType_LP:
	par.pack(buf);
	buf.pack(nl_file_content);
	buf.pack(ipopt_file_content);
	break;
    default:
	abort();
    }
}

/****************************************************************************/

void
BM_lp::unpack_module_data(BCP_buffer& buf)
{
    par.unpack(buf);
    buf.unpack(nl_file_content);
    buf.unpack(ipopt_file_content);

    char* argv_[3];
    char** argv = argv_;
    argv[0] = NULL;
    argv[1] = strdup("dont_even_try_to_open_it.nl");
    argv[2] = NULL;
    std::string ipopt_content(ipopt_file_content.c_str());
    std::string nl_content(nl_file_content.c_str());
    nlp.readAmplNlFile(argv, &ipopt_content, &nl_content, 
                       new Bonmin::IpoptSolver);
    free(argv[1]);

    nlp.extractInterfaceParams();

    /* synchronize bonmin & BCP parameters */
    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
    double bm_intTol;
    double bm_cutoffIncr; // could be negative
    options->GetNumericValue("integer_tolerance",bm_intTol,"bonmin.");
    options->GetNumericValue("cutoff_incr",bm_cutoffIncr,"bonmin.");

    BCP_lp_prob* bcp_lp = getLpProblemPointer();
    const double bcp_intTol = bcp_lp->par.entry(BCP_lp_par::IntegerTolerance);
    const double bcp_cutoffIncr = bcp_lp->par.entry(BCP_lp_par::Granularity);

    if (fabs(bm_intTol - bcp_intTol) > 1e-10) {
	printf("WARNING!\n");
	printf("   The integrality tolerance parameters are different for\n");
	printf("   BCP (%f) and bonmin (%f). They should be identical.\n",
	       bcp_intTol, bm_intTol);
	printf("   For now both will be set to that of bonmin.\n");
    }
    if (fabs(bm_cutoffIncr - bcp_cutoffIncr) > 1e-10) {
	printf("WARNING!\n");
	printf("   The granularity (cutoff increment) parameters are different\n");
	printf("   BCP (%f) and bonmin (%f). They should be identical.\n",
	       bcp_cutoffIncr, bm_cutoffIncr);
	printf("   For now both will be set to that of bonmin.\n");
    }
    bcp_lp->par.set_entry(BCP_lp_par::IntegerTolerance, bm_intTol);
    bcp_lp->par.set_entry(BCP_lp_par::Granularity, bm_cutoffIncr);

    /* If pure BB is selected then a number of BCP parameters are changed */
    if (par.entry(BM_par::PureBranchAndBound)) {
	/* disable strong branching */
	bcp_lp->par.set_entry(BCP_lp_par::MaxPresolveIter, -1);
	/* disable a bunch of printing, all of which are meaningless, since the
	   LP relaxation is meaningless */
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_LpSolutionValue, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_FinalRelaxedSolution, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_RelaxedSolution, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_ReportLocalCutPoolSize, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_ReportLocalVarPoolSize, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_GeneratedCutCount, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_GeneratedVarCount, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_IterationCount, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_RowEffectivenessCount, false);
	//  bcp_lp->par.set_entry(BCP_lp_par::LpVerb_FathomInfo, false);
    }

    /* SOS constraints */
    if (par.entry(BM_par::BranchOnSos)) {
	sos.setFrom(nlp.model()->sosConstraints());
    }

    /* just to be on the safe side... always allocate */
    primal_solution_ = new double[nlp.getNumCols()];
    branching_priority_ = new double[nlp.getNumCols()];

    /* Priorities... Let's see whether the user has set them up. */
    const int* prio = nlp.getPriorities();
    const int numCols = nlp.getNumCols();
    int i;
    if (prio == NULL) {
	// No. Just sett all of them uniformly to 1.0
	for (i = 0; i < numCols; ++i)
	    {
		branching_priority_[i] = 1.0;
	    }
    }
    else {
	/* Yes! the lower the number the earlier we want to branch on it.
	   Let's see how many different numbers there are */
	std::set<int> numbers;
	for (i = 0; i < numCols; ++i)
	    {
		numbers.insert(prio[i]);
	    }
	if (numbers.size() > 15) {
	    /* Sigh... too many... We just have to trust the user's knowledge */
	    if (par.entry(BM_par::CombinedDistanceAndPriority) == true) {
		printf("WARNING!\n");
		printf("   There are too many (>15) different branching priorities\n");
		printf("   defined. Switching off CombinedDistanceAndPriority and\n");
		printf("   will branch on considering only the specified priority\n");
		printf("   order.\n");
		par.set_entry(BM_par::CombinedDistanceAndPriority, false);
	    }
	}
	if (par.entry(BM_par::CombinedDistanceAndPriority) == true) {
	    /* vars with the lowest prio will have their branching_priority_
	       set to 1.0, the next set to 10.0, etc. */
	    const double powerup[15] = {1e14, 1e13, 1e12, 1e11, 1e10, 1e9, 1e8, 1e7,
					1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0};
	    const double powerdo[15] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7,
					1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14};
	    const double* powers =
		par.entry(BM_par::LowPriorityImportant) ? powerup : powerdo;
	    for (i = 0; i < numCols; ++i)
		{
		    const int pos = coinDistance(numbers.begin(), numbers.find(prio[i]));
		    assert (pos >= 0 && pos < 15);
		    branching_priority_[i] = powers[pos];
		}
	}
	else {
	    if (par.entry(BM_par::LowPriorityImportant)) {
		const double maxprio = *numbers.rbegin();
		for (i = 0; i < numCols; ++i)
		    {
			branching_priority_[i] = maxprio - prio[i];
		    }
	    }
	    else {
		for (i = 0; i < numCols; ++i)
		    {
			branching_priority_[i] = prio[i];
		    }
	    }
	}
    }
}

//#############################################################################

void
BM_tm::pack_user_data(const BCP_user_data* ud, BCP_buffer& buf)
{
    const BM_node* data = dynamic_cast<const BM_node*>(ud);
    data->pack(buf);
}

/*---------------------------------------------------------------------------*/

BCP_user_data*
BM_tm::unpack_user_data(BCP_buffer& buf)
{
    return new BM_node(buf);
}

/*****************************************************************************/

void
BM_lp::pack_user_data(const BCP_user_data* ud, BCP_buffer& buf)
{
    const BM_node* data = dynamic_cast<const BM_node*>(ud);
    data->pack(buf);
}

/*---------------------------------------------------------------------------*/

BCP_user_data*
BM_lp::unpack_user_data(BCP_buffer& buf)
{
    BM_node* data = new BM_node(buf);
    numNlpFailed_ = data->numNlpFailed_;
    return data;
}

//#############################################################################
