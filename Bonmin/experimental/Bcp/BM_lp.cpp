// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "OsiClpSolverInterface.hpp"
#include "BM.hpp"
#include "BCP_message_mpi.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_lp.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

// The following is included for "min"
#include "CoinFinite.hpp"

//#############################################################################

BM_lp::BM_lp() :
    BCP_lp_user(),
    babSolver_(3),
    bonmin_(),
    ws_(NULL),
    chooseVar_(NULL),
    freeToBusyRatio_(0.0),
    in_strong(0)
{
      setOsiBabSolver(&babSolver_);
}

/****************************************************************************/

BM_lp::~BM_lp()
{
    delete chooseVar_;
    delete ws_;
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
    BM_node* data = dynamic_cast<BM_node*>(get_user_data());
    numNlpFailed_ = data->numNlpFailed_;

    int i;
    // First copy the bounds into nlp. That way all the branching decisions
    // will be transferred over.
    OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    Bonmin::OsiTMINLPInterface& nlp = *bonmin_.nonlinearSolver();
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
    if (bonmin_.getAlgorithm() > 0) {
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
    if (bonmin_.getAlgorithm() == 0) {
	// if pure B&B
	return test_feasibility_BB(vars);
    } else {
	return test_feasibility_hybrid(lp_result, vars, cuts);
    }
    return NULL; // fake return to quiet gcc
}

/****************************************************************************/

BCP_solution*
BM_lp::test_feasibility_BB(const BCP_vec<BCP_var*>& vars)
{
  Bonmin::OsiTMINLPInterface& nlp = *bonmin_.nonlinearSolver();
    /* PURE_BB TODO:
       Do whatever must be done to:
       1) compute a lower bound using NLP and save it in "lower_bound_"
       2) save the solution of the NLP in "primal_solution_"
       3) if the optimum (well, local opt anyway...) happens to be
       feasible (or along the NLP optimization we have blundered into
       a feas sol) then create "sol"
    */
    BM_solution* sol = NULL;

    char prefix[100];
#ifdef COIN_HAS_MPI
    sprintf(prefix, "%i", getLpProblemPointer()->get_process_id());
#else
    prefix[0] = 0;
#endif
    try {
    switch (par.entry(BM_par::WarmStartStrategy)) {
    case WarmStartNone:
	nlp.initialSolve();
	break;
    case WarmStartFromRoot:
	nlp.setWarmStart(ws_);
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
    }
    catch(Bonmin::TNLPSolver::UnsolvedError &E) {
      E.writeDiffFiles(prefix);
      E.printError(std::cerr);
      //There has been a failure to solve a problem with Ipopt.  And
      //we will output file with information on what has been changed
      //in the problem to make it fail.
      //Now depending on what algorithm has been called (B-BB or
      //other) the failed problem may be at different place.
      //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
    }
    catch(Bonmin::OsiTMINLPInterface::SimpleError &E) {
      std::cerr<<E.className()<<"::"<<E.methodName()
	       <<std::endl
	       <<E.message()<<std::endl;
    }
    catch(CoinError &E) {
      std::cerr<<E.className()<<"::"<<E.methodName()
	       <<std::endl
	       <<E.message()<<std::endl;
    }
    catch (Ipopt::OPTION_INVALID &E) {
      std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
    }
    catch(...) {
      std::cerr<<" unrecognized exception"<<std::endl;
      throw;
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

	for (i = 0; i < numCols; ++i) {
	    if (vars[i]->var_type() == BCP_ContinuousVar)
		continue;
	    const double psol = primal_solution_[i];
	    const double frac = psol - floor(psol);
	    const double dist = min(frac, 1.0-frac);
	    if (dist >= integerTolerance_)
		break;
	}
	if (i == numCols) {
	    /* yipee! a feasible solution! */
	    sol = new BM_solution;
	    //Just copy the solution stored in solver to sol
	    double ptol = 1e-8;
	    // FIXME: I really should use the next line...
	    // nlp.getDblParam(OsiPrimalTolerance, ptol);
	    for (i = 0 ; i < numCols ; i++) {
		if (fabs(primal_solution_[i]) > ptol)
		    sol->add_entry(i, primal_solution_[i]); 
	    }
	    sol->setObjective(lower_bound_);
	}
	if (lower_bound_ > upper_bound()-get_param(BCP_lp_par::Granularity) &&
	    par.entry(BM_par::PrintBranchingInfo)) {
	    printf("BM_lp: At node %i : will fathom because of high lower bound\n",
		   current_index());
	}
    }
    else if (nlp.isProvenPrimalInfeasible()) {
	// prune it!
	// FIXME: if nonconvex, restart from a different place...
	lower_bound_ = 1e200;
	numNlpFailed_ = 0;
	if (par.entry(BM_par::PrintBranchingInfo)) {
	    printf("BM_lp: At node %i : will fathom because of infeasibility\n",
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
    return sol;
}

/****************************************************************************/

BCP_solution*
BM_lp::test_feasibility_hybrid(const BCP_lp_result& lp_result,
			       const BCP_vec<BCP_var*>& vars,
			       const BCP_vec<BCP_cut*>& cuts)
{
    /* First test that the integrality requirements are met. */
    BCP_solution_generic* gsol = test_full(lp_result, vars, integerTolerance_);
    if (! gsol) {}
	
    /* TODO:
       don't test feasibility in every node
    */
    

    OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    
    // Don't test feasibility if the node LP was infeasible
    if (osi->isProvenPrimalInfeasible() ) {   
	return NULL;
    }

    // The babSolver info used is the one containted in osi
    OsiBabSolver * babSolver =
	dynamic_cast<OsiBabSolver *> (osi->getAuxiliaryInfo());
    //assert(babSolver == &babSolver_);
    babSolver->setSolver(*bonmin_.nonlinearSolver()); 
    //Last cut generator is used to check feasibility
    bonmin_.cutGenerators().back().cgl->generateCuts(*osi, cuts_);
    const int numvar = vars.size();
    double* solverSol = new double[numvar];
    double objValue = 1e200;
    
    BM_solution* sol = NULL;
    if (babSolver->solution(objValue, solverSol,numvar)) {
	sol = new BM_solution;
	// Just copy the solution stored in solver to sol
	for (int i = 0 ; i < numvar ; i++) {
	    if (solverSol[i] > lp_result.primalTolerance())
		sol->add_entry(i, solverSol[i]); 
	}
	sol->setObjective(objValue);
    }
    delete[] solverSol;
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
    if (bonmin_.getAlgorithm() > 0) { /* if not pure B&B */
	/* TODO:
	   don't invoke all of them, only the *good* ones. figure out
	   some measurement of how good a generator is.
	*/
	
	OsiSolverInterface& si = *getLpProblemPointer()->lp_solver;
	double rand;
  Bonmin::BabSetupBase::CuttingMethods::const_iterator end = bonmin_.cutGenerators().end();
  end--;//have to go back one element because last cut generator checks feasibility
    
    /* FIXME: fill out the tree info! */
    CglTreeInfo info;
    for(Bonmin::BabSetupBase::CuttingMethods::const_iterator i = bonmin_.cutGenerators().begin() ;
        i != end ; i++){
    rand = 1 - CoinDrand48();
    if(i->frequency > rand){
      i->cgl->generateCuts(si, cuts_, info);
    }
  }
  	
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
    cuts_.dumpCuts();
}

/****************************************************************************/

double
BM_lp::compute_lower_bound(const double old_lower_bound,
                           const BCP_lp_result& lpres,
                           const BCP_vec<BCP_var*>& vars,
                           const BCP_vec<BCP_cut*>& cuts)
{
    double bd;
    if (bonmin_.getAlgorithm() == 0) {
	/* if pure B&B */
	bd = lower_bound_;
    } else {
	bd = CoinMax(babSolver_.mipBound(), lpres.objval());
    }
    return CoinMax(bd, old_lower_bound);
}

/*---------------------------------------------------------------------------*/
