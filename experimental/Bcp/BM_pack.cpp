// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation

#include "BM.hpp"
#include "BCP_lp.hpp"

#include "BonAmplSetup.hpp"

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

    free(argv[1]);

    /* synchronize bonmin & BCP parameters */
    Ipopt::SmartPtr<Ipopt::OptionsList> options = bonmin_.options();

    /** update getting options directly from setup and store them in vars
	local to the BM_lp object */
    integerTolerance_ = bonmin_.getDoubleParameter(BabSetupBase::IntTol);
    // cutOffDecrement_ could be negative
    double cutOffDecrement =
      bonmin_.getDoubleParameter(BabSetupBase::CutoffDecr);

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
    if (fabs(bcp_cutoffIncr - cutOffDecrement) > 1e-10) {
	printf("WARNING!\n");
	printf("   The granularity (cutoff increment) parameters are different\n");
	printf("   BCP (%f) and bonmin (%f). They should be identical.\n",
	       bcp_cutoffIncr, cutOffDecrement);
	printf("   For now both will be set to that of bonmin.\n");
    }
    bcp_lp->par.set_entry(BCP_lp_par::IntegerTolerance, integerTolerance_);
    bcp_lp->par.set_entry(BCP_lp_par::Granularity, cutOffDecrement);

    /* If pure BB is selected then a number of BCP parameters are changed */
    if (bonmin_.getAlgorithm() == 0 /* pure B&B */) {
	/* disable strong branching */
	bcp_lp->par.set_entry(BCP_lp_par::MaxPresolveIter, -1);
	/* disable a bunch of printing, all of which are meaningless, since the
	   LP relaxation is meaningless */
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_ReportLocalCutPoolSize, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_ReportLocalVarPoolSize, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_GeneratedCutCount, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_GeneratedVarCount, false);
	bcp_lp->par.set_entry(BCP_lp_par::LpVerb_RowEffectivenessCount, false);
    }

    /* extract the sos constraints */
    Bonmin::OsiTMINLPInterface& nlp = *bonmin_.nonlinearSolver();
    const Bonmin::TMINLP::SosInfo * sos = nlp.model()->sosConstraints();
    
    int i;
    const int numCols = nlp.getNumCols();
    const double* clb = nlp.getColLower();
    const double* cub = nlp.getColUpper();
    const int* cPrio = nlp.getPriorities();

    /* Find first the integer variables and then the SOS constraints */
    int nObj = 0;
    OsiObject** osiObj = new OsiObject*[numCols + sos->num];
    for (i = 0; i < numCols; ++i) {
	if (nlp.isInteger(i)) {
	    osiObj[nObj] = new OsiSimpleInteger(i, clb[i], cub[i]);
	    if (cPrio) {
	      osiObj[nObj]->setPriority(cPrio[i]);
	    }
	    ++nObj;
	}
    }
#if ! defined(BM_DISREGARD_SOS)
    const int* starts = sos->starts;
    for (i = 0; i < sos->num; ++i) {
	OsiSOS* so = new OsiSOS(NULL, /* FIXME: why does the constr need */
				starts[i+1] - starts[i],
				sos->indices + starts[i],
				sos->weights + starts[i],
				sos->types[i]);
	so->setPriority(sos->priorities ? sos->priorities[i] : 1);
	osiObj[nObj++] = so;
    }
#endif
    nlp.addObjects(nObj, osiObj);

    objNum_ = nObj;
    /* Allocate the storage arrays */
    infInd_ = new int[objNum_];
    infUseful_ = new double[objNum_];
    feasInd_ = new int[objNum_];
    feasUseful_ = new double[objNum_];
    sbResult_ = new BM_SB_result[objNum_];

    /* Sort the objects based on their priority */
    int* prio = new int[objNum_];
    objInd_ = new int[objNum_];
    for (i = 0; i < objNum_; ++i) {
      sbResult_[i].colInd = osiObj[i]->columnNumber();
      objInd_[i] = i;
      prio[i] = osiObj[i]->priority();
    }
    CoinSort_2(prio, prio+objNum_, objInd_);
    delete[] prio;
    
    for (i = 0; i < nObj; ++i) {
	delete osiObj[i];
    }
    delete[] osiObj;
}

//#############################################################################

void
BM_pack::pack_user_data(const BCP_user_data* ud, BCP_buffer& buf)
{
    const BM_node* data = dynamic_cast<const BM_node*>(ud);
    data->pack(buf);
    BM_tm* tm = dynamic_cast<BM_tm*>(user_class);
    if (tm) {
      tm->pack_pseudo_costs(buf);
    }
}

/*---------------------------------------------------------------------------*/

BCP_user_data*
BM_pack::unpack_user_data(BCP_buffer& buf)
{
    BM_node* node = new BM_node(buf);
    BM_lp* lp = dynamic_cast<BM_lp*>(user_class);
    if (lp) {
      lp->unpack_pseudo_costs(buf);
    }
    return node;
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
