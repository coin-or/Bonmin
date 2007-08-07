// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation

#include "BM.hpp"
#include "BCP_lp.hpp"

#include "BonAmplSetup.hpp"

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
    case BCP_ProcessType_TS:
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

    /* PIERRE create the setup for algorithm, last argument indicates that
      continuous relaxation should not be created.*/
    bonmin_.initialize(argv, ipopt_content, nl_content, false);
    bonmin_.nonlinearSolver()->setExposeWarmStart(true);
    babSolver_.setSolver(bonmin_.nonlinearSolver());

    free(argv[1]);
    if (! get_param(BCP_lp_par::MessagePassingIsSerial) &&
	bonmin_.getAlgorithm() == 0 /* pure B&B */ &&
	par.entry(BM_par::WarmStartStrategy) == WarmStartFromParent) {
	printf("\
BM: WarmStartFromParent is not supported for pure B&B in parallel env.\n");
	printf("\
BM: Switching to WarmStartFromRoot.\n");
	par.set_entry(BM_par::WarmStartStrategy, WarmStartFromRoot);
    }

    /* synchronize bonmin & BCP parameters */
    Ipopt::SmartPtr<Ipopt::OptionsList> options = bonmin_.options();

    /** update getting options directly from setup and store them in vars
	local to the BM_lp object */
    integerTolerance_ = bonmin_.getDoubleParameter(BabSetupBase::IntTol);
    // cutOffDecrement_ could be negative
    cutOffDecrement_ = bonmin_.getDoubleParameter(BabSetupBase::CutoffDecr);

    BCP_lp_prob* bcp_lp = getLpProblemPointer();
    const double bcp_intTol = bcp_lp->par.entry(BCP_lp_par::IntegerTolerance);
    const double bcp_cutoffIncr = bcp_lp->par.entry(BCP_lp_par::Granularity);

    if (fabs(integerTolerance_ - bcp_intTol) > 1e-10) {
	printf("WARNING!\n");
	printf("   The integrality tolerance parameters are different for\n");
	printf("   BCP (%f) and bonmin (%f). They should be identical.\n",
	       bcp_intTol, integerTolerance_);
	printf("   For now both will be set to that of bonmin.\n");
    }
    if (fabs(cutOffDecrement_ - cutOffDecrement_) > 1e-10) {
	printf("WARNING!\n");
	printf("   The granularity (cutoff increment) parameters are different\n");
	printf("   BCP (%f) and bonmin (%f). They should be identical.\n",
	       bcp_cutoffIncr, cutOffDecrement_);
	printf("   For now both will be set to that of bonmin.\n");
    }
    bcp_lp->par.set_entry(BCP_lp_par::IntegerTolerance, integerTolerance_);
    bcp_lp->par.set_entry(BCP_lp_par::Granularity, cutOffDecrement_);

    /* Store a few more options in vars local to the BM_lp object */

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
    if (bonmin_.getAlgorithm() == 0 /* pure B&B */) {
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
#if 0 // Pierre cut generators have already been initialized
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
#else
  //Check if OA decomposition is in cut generators and remove it
    for(BabSetupBase::CuttingMethods::iterator i = bonmin_.cutGenerators().begin() ; 
        i != bonmin_.cutGenerators().end() ; i++){
      OACutGenerator2 * oaDec = dynamic_cast<OACutGenerator2 *> (i->cgl);
      if(oaDec)//Disable it
      {
        i->frequency = 0;
      }
    }
#endif
    }

    /* extract the sos constraints */
Bonmin::OsiTMINLPInterface& nlp = *bonmin_.nonlinearSolver();
const Bonmin::TMINLP::SosInfo * sos = nlp.model()->sosConstraints();
    
    int i;
    const int numCols = nlp.getNumCols();
    const double* clb = nlp.getColLower();
    const double* cub = nlp.getColUpper();

    /* Find first the integer variables and then the SOS constraints */
    int nObj = 0;
    OsiObject** osiObj = new OsiObject*[numCols + sos->num];
    for (i = 0; i < numCols; ++i) {
	if (nlp.isInteger(i)) {
	    osiObj[nObj++] = new OsiSimpleInteger(i, clb[i], cub[i]);
	}
    }
    const int* starts = sos->starts;
    for (i = 0; i < sos->num; ++i) {
	OsiSOS* so = new OsiSOS(NULL, /* FIXME: why does the constr need */
				starts[i+1] - starts[i],
				sos->indices + starts[i],
				sos->weights + starts[i],
				sos->types[i]);
	// FIXME: this should go when SOS object can get a priority
	so->setPriority(1);
	osiObj[nObj++] = so;
	
    }
    nlp.addObjects(nObj, osiObj);
    for (i = 0; i < nObj; ++i) {
	delete osiObj[i];
    }
    delete[] osiObj;

    /* just to be on the safe side... always allocate */
    primal_solution_ = new double[nlp.getNumCols()];

    /* solve the initial nlp to get warmstart info in the root */
    nlp.initialSolve();
    ws_ = nlp.getWarmStart();
    if (get_param(BCP_lp_par::MessagePassingIsSerial) &&
	par.entry(BM_par::WarmStartStrategy) == WarmStartFromParent) {
	warmStart[0] = ws_;
	ws_ = NULL;
    }
}

//#############################################################################

void
BM_pack::pack_user_data(const BCP_user_data* ud, BCP_buffer& buf)
{
    const BM_node* data = dynamic_cast<const BM_node*>(ud);
    data->pack(buf);
}

/*---------------------------------------------------------------------------*/

BCP_user_data*
BM_pack::unpack_user_data(BCP_buffer& buf)
{
    return new BM_node(buf);
}

/*---------------------------------------------------------------------------*/

void
BM_pack::pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf)
{
    const BB_cut* bb_cut = dynamic_cast<const BB_cut*>(cut);
    if (!bb_cut)
	throw BCP_fatal_error("pack_cut_algo() : unknown cut type!\n");
    bb_cut->pack(buf);
}

/*---------------------------------------------------------------------------*/
BCP_cut_algo* BM_pack::unpack_cut_algo(BCP_buffer& buf)
{
    return new BB_cut(buf);
}
