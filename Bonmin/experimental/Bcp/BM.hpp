// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006
#ifndef _BM_H
#define _BM_H

#include "BCP_USER.hpp"
#include "BCP_parameters.hpp"
#include "BCP_tm_user.hpp"
#include "BCP_lp_user.hpp"

#include "BB_cut.hpp"

#include "BonIpoptWarmStart.hpp"

#define BM_DISREGARD_SOS

//#############################################################################

class BM_node : public BCP_user_data {
public:
    /** A counter for how many times in a row did the NLP code fail. When the
	NLP fails we branch; hopefully it'll be OK in the children. If it
	fails too many times in a row then we fathom the node: it's hopelessly
	difficult. */
    int numNlpFailed_;
public:
    BM_node() : numNlpFailed_(0) {}
    BM_node(BCP_buffer& buf) : numNlpFailed_(0) {
	buf.unpack(numNlpFailed_);
    }
    ~BM_node() {}

    inline void pack(BCP_buffer& buf) const {
	buf.pack(numNlpFailed_);
    }
};

//#############################################################################
    
enum BM_message {
    BM_StrongBranchRequest,
    BM_StrongBranchResult,
    BM_PseudoCostUpdate
};

enum BM_BoundChange {
    BM_Var_DownBranch,
    BM_Var_UpBranch
};

//#############################################################################
    
class BM_par {
public:
    enum chr_params {
        //
	DisregardPriorities,
	PrintBranchingInfo,
        end_of_chr_params
    };
    enum int_params {
        //
        UsePseudoCosts,
        DecreasingSortInSetupList,
        PreferHighCombinationInBranching,
        NumNlpFailureMax,

	// twice the number of candidates if all candidates have 2 children.
	// We want to do SB on at least this many (if there are this many)
	SBNumBranchesInRoot,
	
	SBNumBranchesInTree,
	// The level where (and below) there are no min number of branches to
	// be considered and SB need not be done (we can use pseudo costs
	// instead) 
	SBMaxLevel,

        end_of_int_params
    };
    enum dbl_params {
        dummy_dbl_param,
        //
        end_of_dbl_params
    };
    enum str_params {
        NL_filename,
        IpoptParamfile,
        //
        end_of_str_params
    };
    enum str_array_params {
        dummy_str_array_param,
        //
        end_of_str_array_params
    };
};

//#############################################################################

class BM_stats {
public:
  BM_stats() :
    numberNodeSolves_(0),
    numberSbSolves_(0),
    numberFixed_(0),
    numberStrongBranching_(0),
    sumStrongBranchingListIndices_(0),
    sumStrongBranchingListPositions_(0.)
  {}

  ~BM_stats();

  inline void incNumberNodeSolves() {
    numberNodeSolves_++;
  }
  inline void incNumberSbSolves(int cnt) {
    numberSbSolves_ += cnt;
  }
  inline void incNumberFixed() {
    numberFixed_++;
  }
  inline void updateStrongBrachingInfo(int chosenIndex, int listLength) {
    numberStrongBranching_++;
    sumStrongBranchingListIndices_ += chosenIndex;
    sumStrongBranchingListPositions_ +=
      (double)(listLength-chosenIndex)/(double)listLength;
  }
private:
  /** Total number of NLP solves as node solves */
  int numberNodeSolves_;
  /** Total number of NLP solves for strong-branching */
  int numberSbSolves_;
  /** Total number of times variables were fixed due to strong branching */
  int numberFixed_;
  /** Total number of times this node did strong branching */
  int numberStrongBranching_;
  /** Sum of all list indices */
  int sumStrongBranchingListIndices_;
  /** Sum of all relative list positions */
  double sumStrongBranchingListPositions_;
};

//#############################################################################

// Data needed to be sent off to do strong branching

struct BM_BranchData {
  // These are the input for doing the branching
  int changeType;
  int objInd;
  int colInd;
  double solval;
  double bd;
  // These are the results of doing the branching
  int status;
  double objval;
  int iter;
  double time;
};

//#############################################################################

class BM_tm : public BCP_tm_user {

public:

    /**@name Private data member */
    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;
    OsiPseudoCosts pseudoCosts_;

public:

    /**@name Constructors and destructors */
    //@{
    /// Default constructor 
    BM_tm() {}

    /// Default destructor
    virtual ~BM_tm() {}
    //@}

    /**@name Packing and unpacking methods */
    //@{
    virtual void pack_module_data(BCP_buffer& buf, BCP_process_t ptype);

    //@}

    /// Pass the core constraints and core variables to bcp
    virtual void initialize_core(BCP_vec<BCP_var_core*>& vars,
				 BCP_vec<BCP_cut_core*>& cuts,
				 BCP_lp_relax*& matrix);

    /** Create the set of extra variables and cuts that should be added to the
	formulation in the root node. Also decide how variable pricing shuld be
	done, that is, if column generation is requested in the
	init_new_phase() method of this class then column
	generation should be performed according to \c pricing_status.
    */
    virtual void
    create_root(BCP_vec<BCP_var*>& added_vars,
		BCP_vec<BCP_cut*>& added_cuts,
		BCP_user_data*& user_data);

    /// Print a feasible solution
    virtual void display_feasible_solution(const BCP_solution* sol);

    /** Process a message that has been sent by another process' user part to
	this process' user part. */
    virtual void
    process_message(BCP_buffer& buf);
  
    void receive_pseudo_cost_update(BCP_buffer& buf);
    void pack_pseudo_costs(BCP_buffer& buf);

    /// Output the final solution
    virtual void display_final_information(const BCP_lp_statistics& lp_stat);

    virtual void init_new_phase(int phase,
				BCP_column_generation& colgen,
				CoinSearchTreeBase*& candidates);

    void readIpopt();

  private:

  /// auxilliary method for handling output for AMPL
  void write_AMPL_solution(const BCP_solution* sol,
			   bool write_file, bool write_screen);

};

//#############################################################################

struct BM_SB_result
{
  /** 0: Not done   1: Only down   2: only up   3: both ways */
  int branchEval; 
  int objInd;
  int colInd;
  int status[2];
  int iter[2];
  double objval[2];
  double varChange[2];
  double time[2];
};

//#############################################################################

#include <OsiAuxInfo.hpp>
#include <OsiCuts.hpp>
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "BonOaFeasChecker.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonEcpCuts.hpp"
#include "BonOACutGenerator2.hpp"

#include "BCP_lp_user.hpp"
#include "BonAmplSetup.hpp"
#include "BonChooseVariable.hpp"

class BM_lp : public BCP_lp_user
{
    /* There's no totalTime_ and nodeTime_. Look at the top of BM.cpp */
    //   double totalTime_;
    //   double nodeTime_;
    int in_strong;

    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;
    BCP_buffer bm_buf;

    /** This contains the setup for running Bonmin in particular nlp
	solver, continuous solver, cut generators,...*/
    Bonmin::BonminAmplSetup bonmin_;

    double integerTolerance_;

    /** A counter for how many times in a row did the NLP code fail. When the
	NLP fails we branch; hopefully it'll be OK in the children. If it
	fails too many times in a row then we fathom the node: it's hopelessly
	difficult. */
    int numNlpFailed_;

    OsiCuts cuts_;

  /** These are the indices of the integral (i.e., things that can be branched
      on) objects in the solver, sorted based on the priority of the
      corresponding objects */
  int* objInd_;
  int objNum_;

  /** Every time when branching decisions are to be made, we create 6 arrays,
      3 for those objects that are infeasible and 3 for those that are
      feasible. infInd_ contains the indices of the objects (into the objects_
      array of objects) that are not feasible, infUseful_ cointains their
      usefulness and infNum_ their number. They are ordered by their priority
      and within that by their usefulness (that depends on their pseudocosts,
      etc.). feasXXX_ contains the same for objects that are feasible, except
      that SOS objects are not listed there, nor variables that are fixed. */
  int*    infInd_;
  double* infUseful_;
  int     infNum_;
  int*    feasInd_;
  double* feasUseful_;
  int     feasNum_;

  /** This is where we keep the results in case of distributed strong
      branching. The length of the array is objNum_ */
  BM_SB_result* sbResult_;
  /** A pointer to the entry that got selected */
  BM_SB_result* bestSbResult_;

  /** The time when we started to process the node */
  double node_start_time;
      
  /** Class for collecting statistics */
  BM_stats bm_stats;

public:
    BM_lp();
    virtual ~BM_lp();

    inline int& numNlpFailed() {
	return (dynamic_cast<BM_node*>(get_user_data()))->numNlpFailed_;
    }

    virtual void
    unpack_module_data(BCP_buffer& buf);

    /** Process a message that has been sent by another process' user part to
	this process' user part. */
    virtual void
    process_message(BCP_buffer& buf);

    virtual OsiSolverInterface *
    initialize_solver_interface();

    virtual void
    load_problem(OsiSolverInterface& osi, BCP_problem_core* core,
		 BCP_var_set& vars, BCP_cut_set& cuts);

    virtual void
    modify_lp_parameters(OsiSolverInterface* lp, bool in_strong_branching);

    virtual BCP_solution*
    test_feasibility(const BCP_lp_result& lp_result,
		     const BCP_vec<BCP_var*>& vars,
		     const BCP_vec<BCP_cut*>& cuts);
    BCP_solution* test_feasibility_BB(const BCP_lp_result& lp_result,
				      const BCP_vec<BCP_var*>& vars);
    BCP_solution* test_feasibility_hybrid(const BCP_lp_result& lp_result,
					  const BCP_vec<BCP_var*>& vars,
					  const BCP_vec<BCP_cut*>& cuts);

    virtual void
    generate_cuts_in_lp(const BCP_lp_result& lpres,
			const BCP_vec<BCP_var*>& vars,
			const BCP_vec<BCP_cut*>& cuts,
			BCP_vec<BCP_cut*>& new_cuts,
			BCP_vec<BCP_row*>& new_rows);
    virtual void
    cuts_to_rows(const BCP_vec<BCP_var*>& vars, // on what to expand
		 BCP_vec<BCP_cut*>& cuts,       // what to expand
		 BCP_vec<BCP_row*>& rows,       // the expanded rows
		 // things that the user can use for lifting cuts if allowed
		 const BCP_lp_result& lpres,
		 BCP_object_origin origin, bool allow_multiple);
    virtual double
    compute_lower_bound(const double old_lower_bound,
			const BCP_lp_result& lpres,
			const BCP_vec<BCP_var*>& vars,
			const BCP_vec<BCP_cut*>& cuts);

    virtual void 
    initialize_new_search_tree_node(const BCP_vec<BCP_var*>& vars,
				    const BCP_vec<BCP_cut*>& cuts,
				    const BCP_vec<BCP_obj_status>& vs,
				    const BCP_vec<BCP_obj_status>& cs,
				    BCP_vec<int>& var_changed_pos,
				    BCP_vec<double>& var_new_bd,
				    BCP_vec<int>& cut_changed_pos,
				    BCP_vec<double>& cut_new_bd);

    virtual BCP_branching_decision
    select_branching_candidates(const BCP_lp_result& lpres,
				const BCP_vec<BCP_var*>& vars,
				const BCP_vec<BCP_cut*>& cuts,
				const BCP_lp_var_pool& local_var_pool,
				const BCP_lp_cut_pool& local_cut_pool,
				BCP_vec<BCP_lp_branching_object*>& cans,
				bool force_branch = false);

    BCP_branching_decision bbBranch(OsiBranchingInformation& brInfo,
				    BCP_vec<BCP_lp_branching_object*>& cands);
    BCP_branching_decision hybridBranch();

    /** Methods invoked from bbBranch() */
    void send_pseudo_cost_update(OsiBranchingInformation& branchInfo);
    void unpack_pseudo_costs(BCP_buffer& buf);
    int sort_objects(OsiBranchingInformation& branchInfo,
		     Bonmin::BonChooseVariable* choose, int& branchNum);
    void clear_SB_results();
    void collect_branch_data(OsiBranchingInformation& branchInfo,
			     OsiSolverInterface* solver,
			     const int branchNum,
			     BM_BranchData* branchData);
    void do_distributed_SB(OsiBranchingInformation& branchInfo,
			   OsiSolverInterface* solver,
			   const CoinWarmStart* cws,
			   const int branchNum,
			   const int* pids, const int pidNum);
    bool isBranchFathomable(int status, double obj);
    int process_SB_results(OsiBranchingInformation& branchInfo,
			   OsiSolverInterface* solver,
			   Bonmin::BonChooseVariable* choose,
			   OsiBranchingObject*& branchObject);
    int try_to_branch(OsiBranchingInformation& branchInfo,
		      OsiSolverInterface* solver,
		      Bonmin::BonChooseVariable* choose,
		      OsiBranchingObject*& branchObject,
		      bool allowVarFix);

    virtual void
    set_user_data_for_children(BCP_presolved_lp_brobj* best, 
			       const int selected);

};

//#############################################################################

#include "BCP_USER.hpp"

class BM_pack : public BCP_user_pack {
public:
    virtual ~BM_pack() {}

    virtual void pack_user_data(const BCP_user_data* ud, BCP_buffer& buf);
    virtual BCP_user_data* unpack_user_data(BCP_buffer& buf);

    virtual void pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf);
    virtual BCP_cut_algo* unpack_cut_algo(BCP_buffer& buf);

};

//#############################################################################

class BM_init : public USER_initialize {

public:

    virtual BCP_tm_user * tm_init(BCP_tm_prob& p,
				  const int argnum,
				  const char * const * arglist);

    virtual BCP_lp_user * lp_init(BCP_lp_prob& p);

    virtual BCP_user_pack * packer_init(BCP_user_class* p);
};

#endif
