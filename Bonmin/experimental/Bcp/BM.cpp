// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006
#include <sstream>

/* To get the cumulative time spent on a processor just use a gawk command
   like this below. Look at the output first; probably the process id needs
   to be prepended to the regexp and the procid may also change the $7 to
   some other word.
   gawk -e 'BEGIN {t=0} /^BCP_lp: Time spent in this node:/ {t+=$7} END {print t}' outfile
*/

#include "CoinHelperFunctions.hpp"
#include "BonAmplInterface.hpp"
#include "BonTMINLP.hpp"
#include "BonIpoptSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include "BCP_parameters.hpp"
#include "BCP_lp.hpp"
#include "BM.hpp"
#include "BB_cut.hpp"
#include "CoinDistance.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_tm.hpp"

using namespace std;

//#############################################################################

int main(int argc, char* argv[])
{
    BM_init user_init;
    return bcp_main(argc, argv, &user_init);
}

//#############################################################################

template <>
void BCP_parameter_set<BM_par>::create_keyword_list() {
    // Create the list of keywords for parameter file reading
    keys.push_back(make_pair(BCP_string("BranchOnSos"),
			     BCP_parameter(BCP_CharPar, BranchOnSos)));
    keys.push_back(make_pair(BCP_string("PureBranchAndBound"),
			     BCP_parameter(BCP_CharPar, PureBranchAndBound)));
    keys.push_back(make_pair(BCP_string("PrintBranchingInfo"),
			     BCP_parameter(BCP_CharPar, PrintBranchingInfo)));
    keys.push_back(make_pair(BCP_string("CombinedDistanceAndPriority"),
			     BCP_parameter(BCP_CharPar, CombinedDistanceAndPriority)));
    keys.push_back(make_pair(BCP_string("LowPriorityImportant"),
			     BCP_parameter(BCP_CharPar, LowPriorityImportant)));
    keys.push_back(make_pair(BCP_string("NumNlpFailureMax"),
			     BCP_parameter(BCP_IntPar, NumNlpFailureMax)));
    keys.push_back(make_pair(BCP_string("NL_filename"),
			     BCP_parameter(BCP_StringPar, NL_filename)));
    keys.push_back(make_pair(BCP_string("IpoptParamfile"),
			     BCP_parameter(BCP_StringPar, IpoptParamfile)));
}

template <>
void BCP_parameter_set<BM_par>::set_default_entries() {
    set_entry(BranchOnSos, true);
    set_entry(PureBranchAndBound, false);
    set_entry(PrintBranchingInfo, true);
    set_entry(CombinedDistanceAndPriority, true);
    set_entry(LowPriorityImportant, true);
    set_entry(NumNlpFailureMax, 5);
    set_entry(NL_filename, "");
    set_entry(IpoptParamfile, "");
}

//#############################################################################

BCP_lp_user *
BM_init::lp_init(BCP_lp_prob& p)
{
    return new BM_lp;
}

/****************************************************************************/

BCP_tm_user *
BM_init::tm_init(BCP_tm_prob& p,
                 const int argnum, const char * const * arglist)
{
    BM_tm* tm = new BM_tm;

    if (argnum == 2) {
	tm->par.read_from_file(arglist[1]);
    } else if (argnum == 1) {
	// work with defaults
    } else {
	tm->par.read_from_arglist(argnum, arglist);
    }

    tm->readIpopt();

    return tm;
}

//#############################################################################

void
BM_tm::readIpopt()
{
    if (par.entry(BM_par::IpoptParamfile).length() != 0) {
	FILE* ipopt = fopen(par.entry(BM_par::IpoptParamfile).c_str(), "r");
	if (!ipopt) {
	    throw BCP_fatal_error("Non-existent IpoptParamfile\n");
	}
	fseek(ipopt, 0, SEEK_END);
	int len = ftell(ipopt) + 1;
	fseek(ipopt, 0, SEEK_SET);
	char* ipopt_content = new char[len+1];
	len = fread(ipopt_content, 1, len, ipopt);
	ipopt_file_content.assign(ipopt_content, len);
	delete[] ipopt_content;
    }

    FILE* nl = fopen(par.entry(BM_par::NL_filename).c_str(), "r");
    if (!nl) {
	throw BCP_fatal_error("Non-existent NL_filename\n");
    }
    fseek(nl, 0, SEEK_END);
    int len = ftell(nl) + 1;
    fseek(nl, 0, SEEK_SET);
    char* nl_content = new char[len+1];
    len = fread(nl_content, 1, len, nl);
    nl_file_content.assign(nl_content, len);
    delete[] nl_content;
}

/****************************************************************************/

void
BM_tm::initialize_core(BCP_vec<BCP_var_core*>& vars,
                       BCP_vec<BCP_cut_core*>& cuts,
                       BCP_lp_relax*& matrix)
{
    Bonmin::AmplInterface nlpSolver; 
    char* argv_[3];
    char** argv = argv_;
    argv[0] = NULL;
    argv[1] = strdup(par.entry(BM_par::NL_filename).c_str());
    argv[2] = NULL;
    nlpSolver.readAmplNlFile(argv, 0, 0, new Bonmin::IpoptSolver);
    free(argv[1]);
  
    nlpSolver.extractInterfaceParams();
  
    OsiClpSolverInterface clp;
    int addObjVar = par.entry(BM_par::PureBranchAndBound) ? 0 : 1;
    nlpSolver.extractLinearRelaxation(clp, addObjVar);
  
    const int numCols = clp.getNumCols();
    const int numRows = clp.getNumRows();

    double* obj = new double[numCols];
    if (par.entry(BM_par::PureBranchAndBound)) {
	CoinFillN(obj, numCols, 0.0);
    }
    else {
	CoinDisjointCopyN(clp.getObjCoefficients(), numCols, obj);
    }

    vars.reserve(numCols);
    const double* clb = clp.getColLower();
    const double* cub = clp.getColUpper();
    for (int i = 0; i < numCols; ++i)
	{
	    BCP_var_t type = BCP_ContinuousVar;
	    if (clp.isBinary(i)) type = BCP_BinaryVar;
	    if (clp.isInteger(i)) type = BCP_IntegerVar;
	    vars.push_back(new BCP_var_core(type, obj[i], clb[i], cub[i]));
	}

    if (par.entry(BM_par::PureBranchAndBound)) {
	// Just fake something into the core matrix. In this case: 0 <= 1
	BCP_vec<double> OBJ(obj, numCols);
	BCP_vec<double> CLB(clb, numCols);
	BCP_vec<double> CUB(cub, numCols);
	BCP_vec<double> RLB(1, 0.0);
	BCP_vec<double> RUB(1, 1.0);
	BCP_vec<int> VB;
	BCP_vec<int> EI;
	BCP_vec<double> EV;
	VB.push_back(0);
	VB.push_back(2);
	EI.push_back(0);         EV.push_back(0.0);
	EI.push_back(numCols-1); EV.push_back(0.0);
	matrix = new BCP_lp_relax(false, VB, EI, EV, OBJ, CLB, CUB, RLB, RUB);
	cuts.push_back(new BCP_cut_core(0.0, 1.0));
    } else {
	cuts.reserve(numRows);
	const double* rlb = clp.getRowLower();
	const double* rub = clp.getRowUpper();
	for (int i = 0; i < numRows; ++i)
	    {
		cuts.push_back(new BCP_cut_core(rlb[i], rub[i]));
	    }

	matrix = new BCP_lp_relax(true /*column major ordered*/);
	matrix->copyOf(*clp.getMatrixByCol(), obj, clb, cub, rlb, rub);
    }
}

/****************************************************************************/
void
BM_tm::create_root(BCP_vec<BCP_var*>& added_vars,
                   BCP_vec<BCP_cut*>& added_cuts,
                   BCP_user_data*& user_data,
                   BCP_pricing_status& pricing_status)
{
    BM_node* data = new BM_node;
    data->numNlpFailed_ = 0;
    user_data = data;
}

/****************************************************************************/

void
BM_tm::display_feasible_solution(const BCP_solution* sol)
{
    const BM_solution* bs = dynamic_cast<const BM_solution*>(sol);
    if (!bs) {
	throw BCP_fatal_error("Trying to pack non-BM_solution.\n");
    }

    /* Parse again the input file so that we have a nice and clean ampl setup */
    Bonmin::AmplInterface nlpSolver;  
    char* argv_[3];
    char** argv = argv_;
    argv[0] = NULL;
    argv[1] = strdup(par.entry(BM_par::NL_filename).c_str());
    argv[2] = NULL;
    nlpSolver.readAmplNlFile(argv, 0, 0, new Bonmin::IpoptSolver);
    free(argv[1]);
    OsiClpSolverInterface clp;
    int addObjVar = par.entry(BM_par::PureBranchAndBound) ? 0 : 1;
    nlpSolver.extractLinearRelaxation(clp, addObjVar);
    const int numCols = clp.getNumCols();

    int i;
  
    /* This will give the vector of core variables we have created in
       BM_tm::initialize_core */
    const BCP_vec<BCP_var_core*>& vars = getTmProblemPointer()->core->vars;

    /* Create a dense vector with the value of each variable. Round the
       integer vars (they were tested to be near enough to integrality) so
       in the printouts we won't have 9.99999991234, etc.) */
    double* dsol = new double[numCols];
    for (i = 0; i < numCols; ++i) {
	dsol[i] = 0.0;
    }
    const int size = bs->_ind.size();
    for (i = 0; i < size; ++i) {
	const int ind = bs->_ind[i];
	const double v = bs->_values[i];
	const BCP_var_t type = vars[ind]->var_type();
	dsol[ind] = (type == BCP_ContinuousVar) ? v : floor(v+0.5);
    }

    /* Display the solution on the screen */
    printf("bonmin: feasible solution found.  Objective value: %f\n",
	   bs->_objective);
    for (i = 0; i < size; ++i) {
	printf("     index: %5i     value: %f\n", bs->_ind[i], dsol[bs->_ind[i]]);
    }
    printf("\n");

    /* create the AMPL solfile */
    nlpSolver.writeAmplSolFile("\nbon-min: Optimal solution", dsol, NULL);
    delete[] dsol;
}

//#############################################################################

BM_lp::BM_lp() :
    BCP_lp_user(),
    babSolver_(3),
    nlp(),
    sos(),
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
    // First copy the bounds into nlp. That way all the branching decisions will
    // be transferred over.
    OsiSolverInterface * osi = getLpProblemPointer()->lp_solver;
    const int numCols = nlp.getNumCols();
    const double* clb = osi->getColLower();
    const double* cub = osi->getColUpper();
    for (int i = 0; i < numCols; ++i) {
	const BCP_var_core* v = dynamic_cast<const BCP_var_core*>(vars[i]);
	if (v) {
	    nlp.setColLower(i, clb[i]);
	    nlp.setColUpper(i, cub[i]);
	    continue;
	}
	const BM_branching_var* bv = dynamic_cast<const BM_branching_var*>(vars[i]);
	if (bv) {
	    const int ind = bv->sos_index;
	    const int split = bv->split;
	    const int len = sos.lengths[ind];
	    const int *indices = sos.indices[ind];
	    if (bv->ub() == 0.0) {
		for (int j = split; j < len; ++j) {
		    nlp.setColLower(indices[j], 0.0);
		    nlp.setColUpper(indices[j], 0.0);
		}
	    } else {
		for (int j = 0; j < split; ++j) {
		    nlp.setColLower(indices[j], 0.0);
		    nlp.setColUpper(indices[j], 0.0);
		}
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
	    lower_bound_ = nlp.getObjValue();
	    CoinDisjointCopyN(nlp.getColSolution(), vars.size(), primal_solution_);
	    numNlpFailed_ = 0;
	    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
	    double intTol;
	    options->GetNumericValue("integer_tolerance",intTol,"bonmin.");

	    const int numCols = nlp.getNumCols();
	    int i;
	    for (i = 0; i < numCols; ++i)
		{
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
		for (i = 0 ; i < numCols ; i++)
		    {
			if (primal_solution_[i] > lp_result.primalTolerance())
			    sol->add_entry(i, primal_solution_[i]); 
		    }
		sol->setObjective(lower_bound_);
	    }
	}
	else if (nlp.isProvenPrimalInfeasible()) {
	    // prune it!
	    lower_bound_ = 1e200;
	    numNlpFailed_ = 0;
	}
	else if (nlp.isAbandoned()) {
	    // nlp failed
	    nlp.forceBranchable();
	    lower_bound_ = nlp.getObjValue();
	    CoinDisjointCopyN(nlp.getColSolution(), vars.size(), primal_solution_);
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
	nlp.retrieve_options()->GetNumericValue("cutoff_incr",
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
	for(int i = 0 ; i < numCuts ; i++)
	    {
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
    for (int i = 0; i < numCuts; ++i)
	{
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

BCP_branching_decision
BM_lp::select_branching_candidates(const BCP_lp_result& lpres,
                                   const BCP_vec<BCP_var*>& vars,
                                   const BCP_vec<BCP_cut*>& cuts,
                                   const BCP_lp_var_pool& local_var_pool,
                                   const BCP_lp_cut_pool& local_cut_pool,
                                   BCP_vec<BCP_lp_branching_object*>& cands)
{
    if (! par.entry(BM_par::PureBranchAndBound)) {
	return BCP_lp_user::select_branching_candidates(lpres, vars, cuts,
							local_var_pool,
							local_cut_pool, cands);
    }

    Ipopt::SmartPtr<Ipopt::OptionsList> options = nlp.retrieve_options();
    double intTol;
    double cutoffIncr;
    options->GetNumericValue("integer_tolerance", intTol, "bonmin.");
    options->GetNumericValue("cutoff_incr", cutoffIncr, "bonmin.");

    if (lower_bound_ >= upper_bound() + cutoffIncr) {
	return BCP_DoNotBranch_Fathomed;
    }

    if (numNlpFailed_ >= par.entry(BM_par::NumNlpFailureMax)) {
	printf("WARNING! Too many (%i) NLP failures in a row. Abandoning node.",
	       numNlpFailed_);
	return BCP_DoNotBranch_Fathomed;
    }

    if (par.entry(BM_par::BranchOnSos) && sos.num > 0) {
	/* FIXME :
	   here we assume that *everything* in an SOS constraint is binary */
	double bestval = 2.0;
	int bestsos = -1;
	int bestsplit = -1;
	int bestprio = 0;
	for (int i = 0; i < sos.num; ++i) {
	    const int ind = sos.priority_order[i];
	    const char type = sos.types[ind];
	    if (type != 1)
		continue;
	    const int prio = sos.priorities[ind];
	    if (bestsos >= 0 && prio != bestprio)
		break;
	    const int len = sos.lengths[ind];
	    const int *indices = sos.indices[ind];
	    double primal = 0.0;
	    for (int j = 0; j < len; ++j) {
		const double prim_i = primal_solution_[indices[j]];
		if (primal + prim_i > 0.5) {
		    if (prim_i < 1-intTol) {
			if (j == 0)
			    ++j;
			/* one branch is [0,j) the other is [j,num) */
			double val = prim_i;
			if (val < bestval) {
			    bestval = val;
			    bestsos = ind;
			    bestsplit = j;
			    bestprio = prio;
			}
		    }
		    break;
		}
		primal += prim_i;
	    }
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
		printf("BM_lp: branching on SOS #%i   split pos: %i\n",
		       bestsos, bestsplit);
	    }
	    return BCP_DoBranch;;
	}
    }

    /* If PURE_BB then we have the NLP optimum in "primal_solution_". 
       Try to find a good branching candidate there. */
    const int numVars = vars.size();
    int besti = -1;
    double bestval = -1e200;
    for (int i = 0; i < numVars; ++i)
	{
	    if (vars[i]->var_type() == BCP_ContinuousVar)
		continue;
	    const double psol = CoinMin(CoinMax(vars[i]->lb(), primal_solution_[i]),
					vars[i]->ub());
	    const double frac = psol - floor(psol);
	    const double dist = std::min(frac, 1.0-frac);

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

BM_lp::BmSosInfo::~BmSosInfo()
{
    if (num > 0) {
	for (int i = 0; i < num; ++i) {
	    delete[] indices[i];
	    delete[] weights[i];
	}
	delete[] indices;
	delete[] weights;
	delete[] types;
	delete[] lengths;
	delete[] priority_order;
    }
}

/*---------------------------------------------------------------------------*/

void BM_lp::BmSosInfo::setFrom(const Bonmin::TMINLP::SosInfo * sos)
{
    if (!sos || sos->num == 0)
	return;

    //we have some sos constraints
    num = sos->num;
    priority_order = new int[num];
    CoinIotaN(priority_order, num, 0);
    priorities = new int[num];
    if (sos->priorities) {
	CoinDisjointCopyN(sos->priorities, num, priorities);
    } else {
	CoinFillN(priorities, num, 0);
    }
    /* FIXME: we may need to go down in order! */
    CoinSort_2(priorities, priorities+num, priority_order);

    lengths = new int[num];
    types = new char[num];
    indices = new int*[num];
    weights = new double*[num];

    for (int i = 0; i < num; ++i) {
	lengths[i] = sos->starts[i+1] - sos->starts[i];
	types[i] = sos->types[i];
	indices[i] = new int[lengths[i]];
	weights[i] = new double[lengths[i]];
	CoinDisjointCopyN(sos->indices+sos->starts[i], lengths[i], indices[i]);
	CoinDisjointCopyN(sos->weights+sos->starts[i], lengths[i], weights[i]);
	CoinSort_2(weights[i], weights[i]+lengths[i], indices[i]);
    }
}

/*---------------------------------------------------------------------------*/

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

//#############################################################################

BCP_MemPool BM_branching_var::memPool(sizeof(BM_branching_var));
