// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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

#ifndef BM_DEBUG_PRINT
#define BM_DEBUG_PRINT 0
#endif

static char prefix[100];

//#############################################################################

BM_lp::BM_lp() :
    BCP_lp_user(),
    in_strong(0),
    bonmin_(),
    objInd_(NULL),
    objNum_(0),
    infInd_(NULL),
    infUseful_(NULL),
    infNum_(0),
    feasInd_(NULL),
    feasUseful_(NULL),
    feasNum_(0),
    sbResult_(NULL),
    bestSbResult_(NULL)
{
}

/****************************************************************************/

BM_lp::~BM_lp()
{
  delete[] objInd_;
  delete[] infInd_;
  delete[] infUseful_;
  delete[] feasInd_;
  delete[] feasUseful_;
  delete[] sbResult_;
}

/****************************************************************************/

OsiSolverInterface *
BM_lp::initialize_solver_interface()
{
#ifdef COIN_HAS_MPI
  sprintf(prefix, "%i", getLpProblemPointer()->get_process_id());
#else
  prefix[0] = 0;
#endif

  OsiSolverInterface* solver = NULL;
  if (bonmin_.getAlgorithm() == 0) {
    // Pure B&B
    solver = bonmin_.nonlinearSolver()->clone();
  } else {
    OsiClpSolverInterface * clp = new OsiClpSolverInterface;
    OsiBabSolver babSolver(3);
    babSolver.setSolver(clp);
    clp->setAuxiliaryInfo(&babSolver);
    clp->messageHandler()->setLogLevel(0);
    setOsiBabSolver(dynamic_cast<OsiBabSolver *>(clp->getAuxiliaryInfo()));
    solver = clp;
  }
  return solver;
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
    node_start_time = CoinWallclockTime();
    BM_node* data = dynamic_cast<BM_node*>(get_user_data());
    numNlpFailed_ = data->numNlpFailed_;

    if (bonmin_.getAlgorithm() != 0) {
      // Not pure BB, so an LP solver will be used. Now we have to...

      // First copy the bounds into nlp. That way all the branching decisions
      // will be transferred over.
      int i;

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

      // copy over the OsiObjects from the nlp solver
      osi->addObjects(nlp.numberObjects(), nlp.objects());
    }

    in_strong = 0;
}

/************************************************************************/
void
BM_lp::load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
		    BCP_var_set& vars, BCP_cut_set& cuts)
{
  if (bonmin_.getAlgorithm() != 0) {
    // We are doing hybrid, so osi is an LP solver. Call the default.
    BCP_lp_user::load_problem(osi, core, vars, cuts);
    return;
  }
  // Otherwise we do B&B and osi is an NLP solver.
  // There is no need to fill it with the data from bonmin_.nonlinearSolver()
  // since osi is a clone() of the master_lp in the BCP_lp object, and the
  // master_lp is a clone() of bonmin_.nonlinearSolver(), and that clone()
  // copies the data as well.
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
    return test_feasibility_BB(lp_result, vars);
  } else {
    return test_feasibility_hybrid(lp_result, vars, cuts);
  }
  return NULL; // fake return to quiet gcc
}

/****************************************************************************/

BCP_solution*
BM_lp::test_feasibility_BB(const BCP_lp_result& lpres,
			   const BCP_vec<BCP_var*>& vars)
{
  // We can just take the primal solution and test whether it satisfies
  // integrality requirements
  BCP_solution_generic* sol = test_full(lpres, vars, integerTolerance_);
  if (sol) {
    sol->set_objective_value(lpres.objval());
#if (BM_DEBUG_PRINT != 0)
    printf("LP %.3f: Solution found. node: %i  depth: %i  value: %f\n",
	   CoinWallclockTime() - start_time(),
	   current_index(), current_level(), lpres.objval());
#endif
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
    babSolver->setSolver(*bonmin_.nonlinearSolver()); 
    //Last cut generator is used to check feasibility
    bonmin_.cutGenerators().back().cgl->generateCuts(*osi, cuts_);
    const int numvar = vars.size();
    double* solverSol = new double[numvar];
    double objValue = 1e200;
    
    BCP_solution_generic* sol = NULL;
    if (babSolver->solution(objValue, solverSol, numvar)) {
        sol = new BCP_solution_generic(false);
	// Just copy the solution stored in solver to sol
	for (int i = 0 ; i < numvar ; i++) {
	    if (solverSol[i] > lp_result.primalTolerance())
		sol->add_entry(vars[i], solverSol[i]); 
	}
	sol->set_objective_value(objValue);
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
	bd = lpres.objval();
    } else {
      /* FIXME: what the heck was this...
	bd = CoinMax(babSolver_.mipBound(), lpres.objval());
      */
    }
    return CoinMax(bd, old_lower_bound);
}

/*---------------------------------------------------------------------------*/
