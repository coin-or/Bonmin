#include "OsiClpSolverInterface.hpp"
#include "BM.hpp"

//#############################################################################

BM_lp::BM_lp() :
    BCP_lp_user(),
    sos(),
    babSolver_(3),
    nlp(),
    feasChecker_(NULL),
    in_strong(0)
{
    babSolver_.setSolver(&nlp);
    setOsiBabSolver(&babSolver_);
}

/****************************************************************************/

BM_lp::~BM_lp()
{
    delete feasChecker_;
    delete[] primal_solution_;
    delete[] branching_priority_;
    for (int i = sos.size() - 1; i >= 0; --i) {
	delete sos[i];
    }
}

/****************************************************************************/

OsiSolverInterface *
BM_lp::initialize_solver_interface()
{
    OsiClpSolverInterface * clp = new OsiClpSolverInterface;
    OsiBabSolver babSolver(3);
    babSolver.setSolver(clp);
    clp->setAuxiliaryInfo(&babSolver);
    clp->messageHandler()->setLogLevel(0);
    setOsiBabSolver(dynamic_cast<OsiBabSolver *>(clp->getAuxiliaryInfo()));
    return clp;
}

/****************************************************************************/

void
BM_lp::initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
				       const BCP_vec<BCP_cut*>& cuts,
				       const BCP_vec<BCP_obj_status>& vs,
				       const BCP_vec<BCP_obj_status>& cs,
				       BCP_vec<int>& var_changed_pos,
				       BCP_vec<double>& var_new_bd,
				       BCP_vec<int>& cut_changed_pos,
				       BCP_vec<double>& cut_new_bd)
{
    int i;
    // First copy the bounds into nlp. That way all the branching decisions
    // will be transferred over.
    for (i = sos.size() - 1; i >= 0; --i) {
	sos[i]->first = 0;
	sos[i]->last = sos[i]->length;
    }
    OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    const int numCols = osi->getNumCols();
    const double* clb = osi->getColLower();
    const double* cub = osi->getColUpper();
    for (i = 0; i < numCols; ++i) {
	const BCP_var_core* v =
	    dynamic_cast<const BCP_var_core*>(vars[i]);
	if (v) {
	    nlp.setColLower(i, clb[i]);
	    nlp.setColUpper(i, cub[i]);
	    continue;
	}
	const BM_branching_var* bv =
	    dynamic_cast<const BM_branching_var*>(vars[i]);
	if (bv) {
	    const int ind = bv->sos_index;
	    const int split = bv->split;
	    const char type = sos[ind]->type;
	    const int *indices = sos[ind]->indices;
	    if (bv->ub() == 0.0) {
		const int last = sos[ind]->last;
		for (int j = split + type; j < last; ++j) {
		    nlp.setColLower(indices[j], 0.0);
		    nlp.setColUpper(indices[j], 0.0);
		}
		sos[ind]->last = split + type;
	    } else {
		const int first = sos[ind]->first;
		for (int j = first; j < split; ++j) {
		    nlp.setColLower(indices[j], 0.0);
		    nlp.setColUpper(indices[j], 0.0);
		}
		sos[ind]->first = split;
	    }
	    continue;
	}
	throw BCP_fatal_error("BM_lp: unrecognized variable type\n");
    }

    in_strong = 0;
}

/************************************************************************/
void
BM_lp::modify_lp_parameters(OsiSolverInterface* lp, bool in_strong_branching)

    // Called each time the node LP is solved

{
    if (in_strong_branching) {
	in_strong = 1;
	//     lp->setIntParam(OsiMaxNumIterationHotStart, 50);
    }
}

/****************************************************************************/

BCP_solution*
BM_lp::test_feasibility(const BCP_lp_result& lp_result,
                        const BCP_vec<BCP_var*>& vars,
                        const BCP_vec<BCP_cut*>& cuts)
{
    BM_solution* sol = NULL;

    if (par.entry(BM_par::PureBranchAndBound)) {
	/* PURE_BB TODO:
	   Do whatever must be done to:
	   1) compute a lower bound using NLP and save it in "lower_bound_"
	   2) save the solution of the NLP in "primal_solution_"
	   3) if the optimum (well, local opt anyway...) happens to be
	   feasible (or along the NLP optimization we have blundered into
	   a feas sol) then create "sol"
	*/
	nlp.initialSolve();
	if (nlp.isProvenOptimal()) {
	    const int numCols = nlp.getNumCols();
	    lower_bound_ = nlp.getObjValue();
	    CoinDisjointCopyN(nlp.getColSolution(), numCols, primal_solution_);
	    numNlpFailed_ = 0;
	    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
	    double intTol;
	    options->GetNumericValue("integer_tolerance",intTol,"bonmin.");

	    int i;
	    for (i = 0; i < numCols; ++i) {
		if (vars[i]->var_type() == BCP_ContinuousVar)
		    continue;
		const double psol = CoinMin(CoinMax(vars[i]->lb(),
						    primal_solution_[i]),
					    vars[i]->ub());
		const double frac = psol - floor(psol);
		const double dist = std::min(frac, 1.0-frac);
		if (dist >= intTol)
		    break;
	    }
	    if (i == numCols) {
		/* yipee! a feasible solution! */
		sol = new BM_solution;
		//Just copy the solution stored in solver to sol
		for (i = 0 ; i < numCols ; i++) {
		    if (primal_solution_[i] > lp_result.primalTolerance())
			sol->add_entry(i, primal_solution_[i]); 
		}
		sol->setObjective(lower_bound_);
	    }
	    if (lower_bound_>upper_bound()-get_param(BCP_lp_par::Granularity)){
		printf("\
BM_lp: At node %i : will fathom because of high lower bound\n",
		       current_index());
	    }
	}
	else if (nlp.isProvenPrimalInfeasible()) {
	    // prune it!
	    lower_bound_ = 1e200;
	    numNlpFailed_ = 0;
	    printf("\
BM_lp: At node %i : will fathom because of infeasibility\n",
		   current_index());
	}
	else if (nlp.isAbandoned()) {
	    printf("\
BM_lp: At node %i : WARNING: nlp is abandoned. Will force branching\n",
		   current_index());
	    // nlp failed
	    nlp.forceBranchable();
	    lower_bound_ = nlp.getObjValue();
	    CoinDisjointCopyN(nlp.getColSolution(), vars.size(),
			      primal_solution_);
	    numNlpFailed_ += 1;
	    // will branch, but save in the user data how many times we have
	    // failed, and if we fail too many times in a row then just abandon
	    // the node.
	}
	else {
	    // complain loudly about a bug
	    throw BCP_fatal_error("Impossible outcome by nlp.initialSolve()\n");
	}
    }
    else {
	/* TODO:
	   don't test feasibility in every node
	*/
	double integerTolerance;
	double cutOffIncrement;
	nlp.retrieve_options()->GetNumericValue("integer_tolerance",
						integerTolerance,"bonmin.");
	/* FIXME: cutoff_decr was called cutoff_incr??? */
	nlp.retrieve_options()->GetNumericValue("cutoff_decr",
						cutOffIncrement,"bonmin.");
	OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    
	// Don't test feasibility if the node LP was infeasible
	if (osi->isProvenPrimalInfeasible() ) {   
	    return sol;
	}

	if (!feasChecker_) {
	    feasChecker_ = new IpCbcOACutGenerator2(&nlp, NULL, NULL, 
						    cutOffIncrement,
						    integerTolerance, 0, 1);
	    feasChecker_->setLocalSearchNodeLimit(0);
	    feasChecker_->setMaxLocalSearch(0);
	    feasChecker_->setMaxLocalSearchPerNode(0);
	}

	// The babSolver info used is the one containted in osi
	OsiBabSolver * babSolver =
	    dynamic_cast<OsiBabSolver *> (osi->getAuxiliaryInfo());
	//assert(babSolver == &babSolver_);
	babSolver->setSolver(nlp); 
	feasChecker_->generateCuts(*osi, cuts_);
	const int numvar = vars.size();
	double* solverSol = new double[numvar];
	double objValue = 1e200;
    
	if (babSolver->solution(objValue, solverSol,numvar)) {
	    sol = new BM_solution;
	    //Just copy the solution stored in solver to sol
	    for (int i = 0 ; i < numvar ; i++) {
		if (solverSol[i] > lp_result.primalTolerance())
		    sol->add_entry(i, solverSol[i]); 
	    }
	    sol->setObjective(objValue);
	}
	delete[] solverSol;
    }
    return sol;
}

/****************************************************************************/

void
BM_lp::generate_cuts_in_lp(const BCP_lp_result& lpres,
                           const BCP_vec<BCP_var*>& vars,
                           const BCP_vec<BCP_cut*>& cuts,
                           BCP_vec<BCP_cut*>& new_cuts,
                           BCP_vec<BCP_row*>& new_rows)
{
    if (! par.entry(BM_par::PureBranchAndBound)) {
	/* TODO:
	   generate cuts with various other Cgl methods (Gomory, etc)
	*/
	// fill cuts
	// eventually fill rows if not in strong branching
	int numCuts = cuts_.sizeRowCuts();
	for(int i = 0 ; i < numCuts ; i++) {
	    const OsiRowCut& cut = cuts_.rowCut(i);
	    new_cuts.push_back(new BB_cut(cut));
	    const CoinPackedVector& row = cut.row();
	    new_rows.push_back(new BCP_row(row, cut.lb(), cut.ub()));
	}
    }
    cuts_.dumpCuts();
}

/****************************************************************************/

void
BM_lp::cuts_to_rows(const BCP_vec<BCP_var*>& vars, // on what to expand
                    BCP_vec<BCP_cut*>& cuts,       // what to expand
                    BCP_vec<BCP_row*>& rows,       // the expanded rows
                    // things that the user can use for lifting cuts if allowed
                    const BCP_lp_result& lpres,
                    BCP_object_origin origin, bool allow_multiple)
{
    const int numCuts = cuts.size();
    for (int i = 0; i < numCuts; ++i) {
	BB_cut* cut = dynamic_cast<BB_cut*>(cuts[i]);
	if (!cut) {
	    throw BCP_fatal_error("Non-\"BB_cut\" in cuts_to_rows!\n");
	}
	const CoinPackedVector& row = cut->row();
	rows.push_back(new BCP_row(row, cuts[i]->lb(), cuts[i]->ub()));
    }
}

/****************************************************************************/

void
BM_lp::vars_to_cols(const BCP_vec<BCP_cut*>& cuts,
		    BCP_vec<BCP_var*>& vars,
		    BCP_vec<BCP_col*>& cols,
		    const BCP_lp_result& lpres,
		    BCP_object_origin origin, bool allow_multiple)
{
    /* All vars better be BM_branching_var, and their columns are empty */
    const int numVars = vars.size();
    for (int i = 0; i < numVars; ++i) {
	BM_branching_var* bv = dynamic_cast<BM_branching_var*>(vars[i]);
	if (!bv) {
	    throw BCP_fatal_error("Non-\"BM_branching_var\" in vars_to_cols!\n");
	}
	cols.push_back(new BCP_col());
    }
}

/****************************************************************************/

double
BM_lp::compute_lower_bound(const double old_lower_bound,
                           const BCP_lp_result& lpres,
                           const BCP_vec<BCP_var*>& vars,
                           const BCP_vec<BCP_cut*>& cuts)
{
    double bd;
    if (par.entry(BM_par::PureBranchAndBound)) {
	bd = lower_bound_;
    } else {
	bd = std::max(babSolver_.mipBound(), lpres.objval());
    }
    return std::max(bd, old_lower_bound);
}

/****************************************************************************/

BM_lp::BmSosInfo::BmSosInfo(const TMINLP::SosInfo * sosinfo, int i)
{
    priority = sosinfo->priorities ? sosinfo->priorities[i] : 0;
    length   = sosinfo->starts[i+1] - sosinfo->starts[i];
    type     = sosinfo->types[i] == 1 ? 0 : 1;
    indices  = new int[length];
    weights  = new double[length];
    first    = 0;
    last     = length;
    CoinDisjointCopyN(sosinfo->indices+sosinfo->starts[i], length, indices);
    CoinDisjointCopyN(sosinfo->weights+sosinfo->starts[i], length, weights);
    CoinSort_2(weights, weights+length, indices);
}
    
/*---------------------------------------------------------------------------*/

void BM_lp::setSosFrom(const TMINLP::SosInfo * sosinfo)
{
    if (!sosinfo || sosinfo->num == 0)
	return;

    //we have some sos constraints
    sos.reserve(sosinfo->num);
    for (int i = 0; i < sosinfo->num; ++i) {
	sos.push_back(new BmSosInfo(sosinfo, i));
    }

    if (sosinfo->priorities) {
	int * priorities = new int[sosinfo->num];
	CoinDisjointCopyN(sosinfo->priorities, sosinfo->num, priorities);
	/* FIXME: we may need to go down in order! */
	if (par.entry(BM_par::SosWithLowPriorityMoreImportant)) {
	    CoinSort_2(priorities, priorities+sosinfo->num, &sos[0],
		       CoinFirstLess_2<int,BM_lp::BmSosInfo*>());
	} else {
	    CoinSort_2(priorities, priorities+sosinfo->num, &sos[0],
		       CoinFirstGreater_2<int,BM_lp::BmSosInfo*>());
	}
	delete[] priorities;
    }
}

/*---------------------------------------------------------------------------*/

void BM_lp::process_message(BCP_buffer& buf)
{
}

/*---------------------------------------------------------------------------*/

#if 0
void BM_lp::BmSosInfo::shuffle()
{
    if (num == 0)
	return;
    int pos = 0;
    while (pos < num) {
	int cnt = 0;
	while (priorities[priority_order[pos]] ==
	       priorities[priority_order[pos+cnt]]) cnt++;
	if (cnt > 1) {
	    std::random_shuffle(priority_order+pos, priority_order+(pos+cnt));
	}
	pos += cnt;
    }
}
#endif
