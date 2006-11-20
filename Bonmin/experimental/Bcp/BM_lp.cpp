#include "OsiClpSolverInterface.hpp"
#include "BM.hpp"
#include "BCP_lp_node.hpp"

//#############################################################################

BM_lp::BM_lp() :
    BCP_lp_user(),
    babSolver_(3),
    nlp(),
    ws(NULL),
    feasChecker_(NULL),
    in_strong(0)
{
    babSolver_.setSolver(&nlp);
    setOsiBabSolver(&babSolver_);
}

/****************************************************************************/

BM_lp::~BM_lp()
{
    delete ws;
    delete feasChecker_;
    delete[] primal_solution_;
    /* FIXME: CHEATING */
    for (std::map<int, CoinWarmStart*>::iterator it = warmStart.begin();
	 it != warmStart.end(); ++it) {
	delete it->second;
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
    OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    nlp.setColLower(osi->getColLower());
    nlp.setColUpper(osi->getColUpper());

    // Carry the changes over to the object lists in nlp
    const int numObj = nlp.numberObjects();
    OsiObject** nlpObj = nlp.objects();
    for (i = 0; i < numObj; ++i) {
	OsiSimpleInteger* io = dynamic_cast<OsiSimpleInteger*>(nlpObj[i]);
	if (io) {
	    io->resetBounds(&nlp);
	} else {
	    // The rest is OsiSOS where we don't need to do anything
	    break;
	}
    }

    // copy over the OsiObjects from the nlp solver if the lp solver is to be
    // used at all (i.e., not pure B&B)
    if (! par.entry(BM_par::PureBranchAndBound)) {
	osi->addObjects(nlp.numberObjects(), nlp.objects());
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
	switch (par.entry(BM_par::WarmStartStrategy)) {
	case WarmStartNone:
	    nlp.initialSolve();
	    break;
	case WarmStartFromRoot:
	    nlp.setWarmStart(ws);
	    nlp.resolve();
	    break;
	case WarmStartFromParent:
	    /* FIXME: CHEAT! this works only in serial mode! */
	    {
		const int ind = current_index();
		const int parentind =
		    getLpProblemPointer()->parent->index;
		std::map<int, CoinWarmStart*>::iterator it =
		    warmStart.find(parentind);
		nlp.setWarmStart(it->second);
		nlp.resolve();
		warmStart[ind] = nlp.getWarmStart();
		bool sibling_seen =  ((ind & 1) == 0) ?
		    warmStart.find(ind-1) != warmStart.end() :
		    warmStart.find(ind+1) != warmStart.end() ;
		if (sibling_seen) {
		    delete it->second;
		    warmStart.erase(it);
		}
	    }
	    break;
	}

	const int numCols = nlp.getNumCols();
	const double* colsol = nlp.getColSolution();
	if (nlp.isProvenOptimal()) {
	    int i;
	    const double* clb = nlp.getColLower();
	    const double* cub = nlp.getColUpper();
	    // Make sure we are within bounds (get rid of rounding problems)
	    for (i = 0; i < numCols; ++i) {
		primal_solution_[i] =
		    CoinMin(CoinMax(clb[i], colsol[i]), cub[i]);
	    }
	    
	    lower_bound_ = nlp.getObjValue();
	    numNlpFailed_ = 0;
	    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
	    double intTol;
	    options->GetNumericValue("integer_tolerance",intTol,"bonmin.");

	    for (i = 0; i < numCols; ++i) {
		if (vars[i]->var_type() == BCP_ContinuousVar)
		    continue;
		const double psol = primal_solution_[i];
		const double frac = psol - floor(psol);
		const double dist = std::min(frac, 1.0-frac);
		if (dist >= intTol)
		    break;
	    }
	    if (i == numCols) {
		/* yipee! a feasible solution! */
		sol = new BM_solution;
		//Just copy the solution stored in solver to sol
		const double ptol = lp_result.primalTolerance();
		for (i = 0 ; i < numCols ; i++) {
		    if (fabs(primal_solution_[i]) > ptol)
			sol->add_entry(i, primal_solution_[i]); 
		}
		sol->setObjective(lower_bound_);
	    }
	    if (lower_bound_ > upper_bound()-get_param(BCP_lp_par::Granularity) &&
		par.entry(BM_par::PrintBranchingInfo)) {
		printf("\
BM_lp: At node %i : will fathom because of high lower bound\n",
		       current_index());
	    }
	}
	else if (nlp.isProvenPrimalInfeasible()) {
	    // prune it!
	    // FIXME: if nonconvex, restart from a different place...
	    lower_bound_ = 1e200;
	    numNlpFailed_ = 0;
	    if (par.entry(BM_par::PrintBranchingInfo)) {
		printf("\
BM_lp: At node %i : will fathom because of infeasibility\n",
		       current_index());
	    }
	}
	else if (nlp.isAbandoned()) {
	    if (nlp.isIterationLimitReached()) {
		printf("\
BM_lp: At node %i : WARNING: nlp reached iter limit. Will force branching\n",
		       current_index());
	    } else {
		printf("\
BM_lp: At node %i : WARNING: nlp is abandoned. Will force branching\n",
		       current_index());
	    }
	    // nlp failed
	    nlp.forceBranchable();
	    lower_bound_ = nlp.getObjValue();
	    CoinDisjointCopyN(colsol, numCols, primal_solution_);
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
    else { // NOT pure B&B
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
	    feasChecker_ = new Bonmin::OACutGenerator2(&nlp, NULL, NULL,
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
	bd = CoinMax(babSolver_.mipBound(), lpres.objval());
    }
    return CoinMax(bd, old_lower_bound);
}

/*---------------------------------------------------------------------------*/

void BM_lp::process_message(BCP_buffer& buf)
{
}

/*---------------------------------------------------------------------------*/
