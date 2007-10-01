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
    
// This needs to be the same as enum VarSelectStra_Enum in
// BonOsiTMINLPInterface.hpp
enum BM_BranchingStrategy {
    BM_MostFractional=0,
    BM_StrongBranching,
    BM_ReliabilityBranching,
    BM_CurvatureEstimator,
    BM_QpStrongBranching,
    BM_LpStrongBranching,
    BM_NlpStrongBranching,
    BM_OsiChooseVariable,
    BM_OsiChooseStrong
};

enum BM_message {
    StrongBranchRequest,
    StrongBranchResult,
    FreeToBusyRatio
};

//#############################################################################
    
class BM_par {
public:
    enum chr_params {
        //
        CombinedDistanceAndPriority,
	PrintBranchingInfo,
	SosWithLowPriorityMoreImportant,
        VarWithLowPriorityMoreImportant,
        end_of_chr_params
    };
    enum int_params {
        //
	BranchingStrategy,
	FullStrongBranch,
        NumNlpFailureMax,
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

class BM_tm : public BCP_tm_user {

public:

    /**@name Private data member */
    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;

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

class BM_lp : public BCP_lp_user
{
    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;

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

    virtual void
    set_user_data_for_children(BCP_presolved_lp_brobj* best, 
			       const int selected);

private:
    /* There's no totalTime_ and nodeTime_. Look at the top of BM.cpp */
    //   double totalTime_;
    //   double nodeTime_;
    int in_strong;

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
