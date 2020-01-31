// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "OsiClpSolverInterface.hpp"

#include "BCP_problem_core.hpp"
#include "BCP_tm.hpp"
#include "BCP_lp.hpp"

// #include "BonIpoptSolver.hpp"
// #ifdef COIN_HAS_FILTERSQP
// # include "BonFilterSolver.hpp"
// #endif
// #include "BonOACutGenerator2.hpp"
// #include "BonEcpCuts.hpp"
// #include "BonOaNlpOptim.hpp"
#include "BonAmplSetup.hpp"

#include "BM.hpp"

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
    
    /* Get the basic options. */
    Bonmin::BonminAmplSetup bonmin;
    bonmin.readOptionsFile(par.entry(BM_par::IpoptParamfile).c_str());
    bonmin.initialize(argv);    
    
    free(argv[1]);

    Bonmin::OsiTMINLPInterface& nlp = *bonmin.nonlinearSolver();
#if defined(BM_DISREGARD_SOS)
    const Bonmin::TMINLP::SosInfo * sos = nlp.model()->sosConstraints();
    if (sos->num > 0) {
      printf("There are SOS constraints... disregarding them...\n");
    }
#endif
    
    OsiSolverInterface& clp  = *bonmin.continuousSolver();
    
    const int numCols = clp.getNumCols();
    const int numRows = clp.getNumRows();

    const double* clb = clp.getColLower();
    const double* cub = clp.getColUpper();

    double* obj = new double[numCols];
    if (bonmin.getAlgorithm() == Bonmin::B_BB /* pure B&B */) {
      std::cout<<"Doing branch and bound"<<std::endl;
      CoinFillN(obj, numCols, 0.0);
      matrix = NULL;
    } else {
      std::cout<<"Doing hybrid"<<std::endl;
      CoinDisjointCopyN(clp.getObjCoefficients(), numCols, obj);
      cuts.reserve(numRows);
      const double* rlb = clp.getRowLower();
      const double* rub = clp.getRowUpper();
      for (int i = 0; i < numRows; ++i)	{
	cuts.push_back(new BCP_cut_core(rlb[i], rub[i]));
      }
      matrix = new BCP_lp_relax(true /*column major ordered*/);
      matrix->copyOf(*clp.getMatrixByCol(), obj, clb, cub, rlb, rub);
    }

    vars.reserve(numCols);
    for (int i = 0; i < numCols; ++i)	{
      BCP_var_t type = BCP_ContinuousVar;
      if (clp.isBinary(i)) type = BCP_BinaryVar;
      if (clp.isInteger(i)) type = BCP_IntegerVar;
      vars.push_back(new BCP_var_core(type, obj[i], clb[i], cub[i]));
    }
    delete[] obj;

    /* Initialize pseudocosts */
    int nObj = 0;
    for (int i = 0; i < numCols; ++i) {
	if (nlp.isInteger(i)) {
	    ++nObj;
	}
    }
#if ! defined(BM_DISREGARD_SOS)
    const Bonmin::TMINLP::SosInfo * sos = nlp.model()->sosConstraints();
    nObj += sos->num;
#endif
    pseudoCosts_.initialize(nObj);
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
  const BCP_solution_generic* bg =
    dynamic_cast<const BCP_solution_generic*>(sol);

  /* Parse again the input file so that we have a nice and clean ampl
     setup */  
  char* argv_[3];
  char** argv = argv_;
  argv[0] = NULL;
  argv[1] = strdup(par.entry(BM_par::NL_filename).c_str());
  argv[2] = NULL;
  Bonmin::BonminAmplSetup bonmin;
  bonmin.initialize(argv);    
  
  Bonmin::OsiTMINLPInterface& nlpSolver = *bonmin.nonlinearSolver();
  free(argv[1]);
  OsiSolverInterface& clp = *bonmin.continuousSolver();

  const int numCols = clp.getNumCols();

  int i;
  
  /* Create a dense vector with the value of each variable. Round the
     integer vars (they were tested to be near enough to integrality) so
     in the printouts we won't have 9.99999991234, etc.) */
  double* dsol = new double[numCols];
  for (i = 0; i < numCols; ++i) {
    dsol[i] = 0.0;
  }
  const BCP_vec<BCP_var*>& vars = bg->_vars;
  const BCP_vec<double>& values = bg->_values;;
  const int size = vars.size();
  for (i = 0; i < size; ++i) {
    const int ind = vars[i]->bcpind();
    const double v = values[i];
    const BCP_var_t type = vars[i]->var_type();
    dsol[ind] = (type == BCP_ContinuousVar) ? v : floor(v+0.5);
  }

  if (write_screen) {
    /* Display the solution on the screen */
    printf("bonmin: feasible solution found.  Objective value: %f\n",
	   bg->objective_value());
    for (i = 0; i < size; ++i) {
      printf("    index: %5i   value: %f\n",
	     vars[i]->bcpind(), dsol[vars[i]->bcpind()]);
    }
    printf("\n");
  }

  if (write_file) {
    /* create the AMPL solfile */
    nlpSolver.model()->finalize_solution(Bonmin::TMINLP::SUCCESS, nlpSolver.getNumCols(), 
                                         dsol, bg->objective_value());
  }
  delete[] dsol;
}

//#############################################################################

void
BM_tm::process_message(BCP_buffer& buf)
{
  int msgtag;
  buf.unpack(msgtag);
  if (msgtag == BM_PseudoCostUpdate) {
    receive_pseudo_cost_update(buf);
  }
}

/****************************************************************************/

void
BM_tm::display_final_information(const BCP_lp_statistics& lp_stat)
{
  bool write_screen = false;
  BCP_tm_prob *p = getTmProblemPointer();
  if (p->param(BCP_tm_par::TmVerb_FinalStatistics)) {
    printf("TM: Running time: %.3f\n", CoinWallclockTime() - p->start_time);
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
  static const std::string compName;
  static inline const char* name() { return "BMSearchTreeCompareBest"; }
  inline bool operator()(const CoinTreeSiblings* x,
			 const CoinTreeSiblings* y) const {
    double allowable_gap = 1e-8;
    const double xq = x->currentNode()->getQuality();
    const double yq = y->currentNode()->getQuality();
    const int xd = x->currentNode()->getDepth();
    const int yd = y->currentNode()->getDepth();
    return ((xq < yq-allowable_gap) ||
	    (xq <= yq+allowable_gap && xd > yd));
  }
};

const std::string
BMSearchTreeCompareBest::compName("BMSearchTreeCompareBest");

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
    printf("Creating candidate list with BMSearchTreeCompareBest\n");
    break;
  case BCP_BreadthFirstSearch:
    candidates = new CoinSearchTree<CoinSearchTreeCompareBreadth>;
    printf("Creating candidate list with CoinSearchTreeCompareBreadth\n");
    break;
  case BCP_DepthFirstSearch:
    candidates = new CoinSearchTree<CoinSearchTreeCompareDepth>;
    printf("Creating candidate list with CoinSearchTreeCompareDepth\n");
    break;
  case BCP_PreferredFirstSearch:
    candidates = new CoinSearchTree<CoinSearchTreeComparePreferred>;
    printf("Creating candidate list with CoinSearchTreeComparePreferred\n");
    break;
  }
}

//#############################################################################

void
BM_tm::receive_pseudo_cost_update(BCP_buffer& buf)
{
  double* upTotalChange = pseudoCosts_.upTotalChange();
  double* downTotalChange = pseudoCosts_.downTotalChange();
  int* upNumber = pseudoCosts_.upNumber();
  int* downNumber = pseudoCosts_.downNumber();
  int objInd;
  int branch;
  double pcChange;
  while (true) {
    buf.unpack(objInd);
    if (objInd == -1) {
      break;
    }
    buf.unpack(branch);
    buf.unpack(pcChange);
    if (branch == 0) {
      downTotalChange[objInd] += pcChange;
      ++downNumber[objInd];
    } else {
      upTotalChange[objInd] += pcChange;
      ++upNumber[objInd];
    }
  }
}

//-----------------------------------------------------------------------------

void
BM_tm::pack_pseudo_costs(BCP_buffer& buf)
{
  buf.pack(pseudoCosts_.upTotalChange(), pseudoCosts_.numberObjects());
  buf.pack(pseudoCosts_.upNumber(), pseudoCosts_.numberObjects());
  buf.pack(pseudoCosts_.downTotalChange(), pseudoCosts_.numberObjects());
  buf.pack(pseudoCosts_.downNumber(), pseudoCosts_.numberObjects());
}
