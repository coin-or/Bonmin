#include "OsiClpSolverInterface.hpp"
// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "BonIpoptSolver.hpp"
#ifdef COIN_HAS_FILTERSQP
# include "BonFilterSolver.hpp"
#endif
#include "BCP_problem_core.hpp"
#include "BM.hpp"


#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

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
    char* argv_[3];
    char** argv = argv_;
    argv[0] = NULL;
    argv[1] = strdup(par.entry(BM_par::NL_filename).c_str());
    argv[2] = NULL;
    
    Bonmin::AmplInterface nlpSolver(argv);
    Bonmin::OACutGenerator2::registerOptions(nlpSolver.regOptions());
    Bonmin::EcpCuts::registerOptions(nlpSolver.regOptions());
    Bonmin::OaNlpOptim::registerOptions(nlpSolver.regOptions());
    
    free(argv[1]);
    nlpSolver.extractInterfaceParams();
  
    OsiClpSolverInterface clp;
    int addObjVar = minlpParams_.algo == 0 /* pure B&B */ ? 0 : 1;
    nlpSolver.extractLinearRelaxation(clp, addObjVar);
  
    const int numCols = clp.getNumCols();
    const int numRows = clp.getNumRows();

    double* obj = new double[numCols];
    if (minlpParams_.algo == 0 /* pure B&B */) {
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

    if (minlpParams_.algo == 0 /* pure B&B */) {
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
    delete[] obj;
}

/****************************************************************************/
void
BM_tm::create_root(BCP_vec<BCP_var*>& added_vars,
                   BCP_vec<BCP_cut*>& added_cuts,
                   BCP_user_data*& user_data)
{
    BM_node* data = new BM_node;
    data->numNlpFailed_ = 0;
    user_data = data;
}

void
BM_tm::write_AMPL_solution(const BCP_solution* sol,
			   bool write_file, bool write_screen)
{
  const BM_solution* bs = dynamic_cast<const BM_solution*>(sol);
  if (!bs) {
    throw BCP_fatal_error("Trying to pack non-BM_solution.\n");
  }
  /* Parse again the input file so that we have a nice and clean ampl
     setup */
  Bonmin::AmplInterface nlpSolver;  
  Bonmin::OACutGenerator2::registerOptions(nlpSolver.regOptions());
  Bonmin::EcpCuts::registerOptions(nlpSolver.regOptions());
  Bonmin::OaNlpOptim::registerOptions(nlpSolver.regOptions());

  char* argv_[3];
  char** argv = argv_;
  argv[0] = NULL;
  argv[1] = strdup(par.entry(BM_par::NL_filename).c_str());
  argv[2] = NULL;
  nlpSolver.readAmplNlFile(argv, 0, 0);
  minlpParams_.extractParams(&nlpSolver);
  free(argv[1]);
  OsiClpSolverInterface clp;
  int addObjVar = minlpParams_.algo == 0 /* pure B&B */ ? 0 : 1;
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

  if (write_screen) {
    /* Display the solution on the screen */
    printf("bonmin: feasible solution found.  Objective value: %f\n",
	   bs->_objective);
    for (i = 0; i < size; ++i) {
      printf("    index: %5i   value: %f\n", bs->_ind[i], dsol[bs->_ind[i]]);
    }
    printf("\n");
  }

  if (write_file) {
    /* create the AMPL solfile */
    nlpSolver.writeAmplSolFile("\nbon-min: Optimal solution", dsol, NULL);
  }
  delete[] dsol;
}

//#############################################################################

const BCP_proc_id*
BM_tm::process_id() const
{
}

void
BM_tm::send_message(const BCP_proc_id* const target, const BCP_buffer& buf)
{
}

void
BM_tm::broadcast_message(const BCP_process_t proc_type, const BCP_buffer& buf)
{
}

void
BM_tm::process_message(BCP_buffer& buf)
{
}

/****************************************************************************/

void
BM_tm::display_final_information(const BCP_lp_statistics& lp_stat)
{
  bool write_screen = false;
  BCP_tm_prob *p = getTmProblemPointer();
  if (p->param(BCP_tm_par::TmVerb_FinalStatistics)) {
    printf("TM: Running time: %.3f\n", CoinCpuTime() - p->start_time);
    printf("TM: search tree size: %i   ( processed %i )   max depth: %i\n",
	   int(p->search_tree.size()), int(p->search_tree.processed()),
	   p->search_tree.maxdepth());
    lp_stat.display();

    if (! p->feas_sol) {
      printf("TM: No feasible solution is found\n");
    } else {
      printf("TM: The best solution found has value %f\n",
	     p->feas_sol->objective_value());
      if (p->param(BCP_tm_par::TmVerb_BestFeasibleSolution)) {
	write_screen = true;
      }
    }
  }
  if (p->feas_sol) {
    write_AMPL_solution(p->feas_sol, true, write_screen);
  }
}

/****************************************************************************/

void
BM_tm::display_feasible_solution(const BCP_solution* sol)
{
  // For now, we want to write also the AMPL solution file...(?)
  // This needs to be changed though
  write_AMPL_solution(sol, true, true);
}

struct BMSearchTreeCompareBest {
  inline bool operator()(const CoinTreeSiblings* x,
			 const CoinTreeSiblings* y) const {
    double allowable_gap = 1e-8;
    return ((x->currentNode()->quality_ <
	     y->currentNode()->quality_-allowable_gap) ||
	    ((x->currentNode()->quality_ <= y->currentNode()->quality_+allowable_gap) &&
	     x->currentNode()->depth_ > y->currentNode()->depth_));
    }
};

void
BM_tm::init_new_phase(int phase,
		      BCP_column_generation& colgen,
		      CoinSearchTreeBase*& candidates)
{
  BCP_tm_prob* p = getTmProblemPointer();
  colgen = BCP_DoNotGenerateColumns_Fathom;
  switch (p->param(BCP_tm_par::TreeSearchStrategy)) {
  case BCP_BestFirstSearch:
    candidates = new CoinSearchTree<BMSearchTreeCompareBest>;
    break;
  case BCP_BreadthFirstSearch:
    candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
    break;
  case BCP_DepthFirstSearch:
    candidates = new CoinSearchTree<CoinSearchTreeCompareDepth>;
    break;
  }
}

