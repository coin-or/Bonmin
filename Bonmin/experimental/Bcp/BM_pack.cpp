// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
//

#include "BonCbc.hpp"
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
    using namespace Bonmin;

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
    nlp_.readAmplNlFile(argv, NULL, &ipopt_content, &nl_content);
    free(argv[1]);

    nlp_.extractInterfaceParams();
    minlpParams_.extractParams(&nlp_);

    if (! get_param(BCP_lp_par::MessagePassingIsSerial) &&
	minlpParams_.algo == 0 /* pure B&B */ &&
	par.entry(BM_par::WarmStartStrategy) == WarmStartFromParent) {
	printf("\
BM: WarmStartFromParent is not supported for pure B&B in parallel env.\n");
	printf("\
BM: Switching to WarmStartFromRoot.\n");
	par.set_entry(BM_par::WarmStartStrategy, WarmStartFromRoot);
    }

    /* synchronize bonmin & BCP parameters */
    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp_.retrieve_options();

    int nlpLogLevel;
    options->GetIntegerValue("nlp_log_level", nlpLogLevel, "bonmin.");
    nlp_.messageHandler()->setLogLevel(nlpLogLevel);

    double bm_intTol;
    double bm_cutoffIncr; // could be negative
    options->GetNumericValue("integer_tolerance",bm_intTol,"bonmin.");
    options->GetNumericValue("cutoff_decr",bm_cutoffIncr,"bonmin.");

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

    /* Store a few options in local variables */

    options->GetNumericValue("integer_tolerance", integerTolerance_,"bonmin.");
    options->GetNumericValue("cutoff_decr", cutOffDecrement_,"bonmin.");

    // Getting the options for the choose variable object
    if (!options->GetEnumValue("varselect_stra",varselect_,"bonmin.")) {
      // For Bcp, we change the default to most-fractional for now
      varselect_ = Bonmin::OsiTMINLPInterface::MOST_FRACTIONAL;
    }
    options->GetIntegerValue("number_ecp_rounds", numEcpRounds_,"bonmin.");
    options->GetIntegerValue("number_strong_branch",numberStrong_,"bonmin.");
    options->GetIntegerValue("number_before_trust", minReliability_,"bonmin.");
    delete chooseVar_;
    chooseVar_ = NULL;


    /* If pure BB is selected then a number of BCP parameters are changed */
    if (minlpParams_.algo == 0 /* pure B&B */) {
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
    } else {
	/* for hybrid: initialize the cut generators */

	/* NOTE:
	   
	   if the localSerchSolver for oaDec_ is NULL, that will force the cut
	   generator to create a clone of the solverinterface in bcp and use
	   that to solve the milp subproblem.

	   If localSearchSolver is set then if it is *not* the same as bcp's
	   solverinterface then bcp's si's data gets copied into
	   localsearchsolver when the milp is solvede.

	   Finally, if localSearchSolver is set and it is the same as bcp's
	   si, then bcp's si will be modified by the oaDec_. A big NO-NO!
	   Fortunately, this will never happen as in each search tree node bcp
	   creates a clone of the master lp.
	*/
	OsiSolverInterface * localSearchSolver = NULL;
	if (minlpParams_.milpSubSolver == 2) {/* try to use cplex */
#ifdef COIN_HAS_CPX
	    localSearchSolver = new OsiCpxSolverInterface;
	    nlpSolver->extractLinearRelaxation(*localSearchSolver);
#else
	    std::cerr << "You have set an option to use CPLEX as the milp\n"
		      << "subsolver in oa decomposition. However, apparently\n"
		      << "CPLEX is not configured to be used in bonmin.\n"
		      << "See the manual for configuring CPLEX\n";
	    throw -1;
#endif
	}
	Bonmin::initializeCutGenerators(minlpParams_, &nlp_,
					miGGen_, probGen_,
					knapsackGen_, mixedGen_,
					oaGen_, ecpGen_,
					oaDec_, localSearchSolver,
					feasCheck_, NULL
					);
    }

    /* extract the sos constraints */
    const Bonmin::TMINLP::SosInfo * sos = nlp_.model()->sosConstraints();
    
    int i;
    const int numCols = nlp_.getNumCols();
    const double* clb = nlp_.getColLower();
    const double* cub = nlp_.getColUpper();

    /* Find first the integer variables and then the SOS constraints */
    int nObj = 0;
    OsiObject** osiObj = new OsiObject*[numCols + sos->num];
    for (i = 0; i < numCols; ++i) {
	if (nlp_.isInteger(i)) {
	    osiObj[nObj++] = new OsiSimpleInteger(i, clb[i], cub[i]);
	}
    }
    const int* starts = sos->starts;
    for (i = 0; i < sos->num; ++i) {
	osiObj[nObj++] = new OsiSOS(NULL, /* FIXME: why does the constr need */
				    starts[i+1] - starts[i],
				    sos->indices + starts[i],
				    sos->weights + starts[i],
				    sos->types[i]);
    }
    nlp_.addObjects(nObj, osiObj);
    for (i = 0; i < nObj; ++i) {
	delete osiObj[i];
    }
    delete[] osiObj;

    /* just to be on the safe side... always allocate */
    primal_solution_ = new double[nlp_.getNumCols()];

    /* solve the initial nlp to get warmstart info in the root */
    nlp_.initialSolve();
    ws_ = nlp_.getWarmStart();
    if (get_param(BCP_lp_par::MessagePassingIsSerial) &&
	par.entry(BM_par::WarmStartStrategy) == WarmStartFromParent) {
	warmStart[0] = ws_;
	ws_ = NULL;
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
