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

#include "BonCbcParam.hpp"

#include "BCP_USER.hpp"
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
    
enum BM_WarmStartStrategy {
    WarmStartNone,
    WarmStartFromRoot,
    WarmStartFromParent
};

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
	WarmStartStrategy,
        end_of_int_params
    };
    enum dbl_params {
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
    Bonmin::BonminCbcParam minlpParams_;

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

    virtual void pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf) {
	BB_pack_cut(cut, buf);
    }
    virtual BCP_cut_algo* unpack_cut_algo(BCP_buffer& buf) {
	return BB_unpack_cut(buf);
    }

    virtual BCP_solution* unpack_feasible_solution(BCP_buffer& buf);

    /// Packing of user data
    virtual void pack_user_data(const BCP_user_data* ud, BCP_buffer& buf);

    /// Unpacking of user_data
    virtual BCP_user_data* unpack_user_data(BCP_buffer& buf);
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

    /** What is the process id of the current process */
    const BCP_proc_id*
    process_id() const;
    /** Send a message to a particular process */
    void
    send_message(const BCP_proc_id* const target, const BCP_buffer& buf);
    /** Broadcast the message to all processes of the given type */
    void
    broadcast_message(const BCP_process_t proc_type, const BCP_buffer& buf);
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
#include "BonAmplInterface.hpp"
#include "BonTMINLP.hpp"

class BM_lp : public BCP_lp_user
{
    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;

    OsiBabSolver babSolver_;
    Bonmin::AmplInterface nlp_;
    Bonmin::BonminCbcParam minlpParams_;

    CoinWarmStart* ws_;
    OsiChooseVariable* chooseVar_;
    int numEcpRounds_;
    int numberStrong_;
    int minReliability_;
    int varselect_;
    double integerTolerance_;
    double cutOffDecrement_;

    /* A couple of cut generators to be used in the hybrid method */
    CglGomory miGGen_;
    CglProbing probGen_;
    CglKnapsackCover knapsackGen_;
    CglMixedIntegerRounding mixedGen_;
    Bonmin::OaNlpOptim oaGen_;
    Bonmin::EcpCuts ecpGen_;
    Bonmin::OACutGenerator2 oaDec_;
    Bonmin::OaFeasibilityChecker feasCheck_;

    /* FIXME: gross cheating. works only for serial mode. Store the warmstart
       informations in the lp process, do not send them over in user data or
       anywhere. MUST be fixed. The map is indexed by the node index. */
    std::map<int, CoinWarmStart*> warmStart;

    double lower_bound_;
    double* primal_solution_;

    /** A counter for how many times in a row did the NLP code fail. When the
	NLP fails we branch; hopefully it'll be OK in the children. If it
	fails too many times in a row then we fathom the node: it's hopelessly
	difficult. */
    int numNlpFailed_;

    OsiCuts cuts_;

    /** The last free-to-busy ratio (among the LP processes) that the TM has
	sent over. It is used to decide whether we want distributed strong
	branching */
    double freeToBusyRatio_;

public:
    BM_lp();
    virtual ~BM_lp();

    inline int& numNlpFailed() {
	return (dynamic_cast<BM_node*>(get_user_data()))->numNlpFailed_;
    }

    virtual void
    unpack_module_data(BCP_buffer& buf);

    virtual void pack_cut_algo(const BCP_cut_algo* cut, BCP_buffer& buf) {
	BB_pack_cut(cut, buf);
    }
    virtual BCP_cut_algo* unpack_cut_algo(BCP_buffer& buf) {
	return BB_unpack_cut(buf);
    }

    virtual void
    pack_feasible_solution(BCP_buffer& buf, const BCP_solution* sol);

    virtual void
    pack_user_data(const BCP_user_data* ud, BCP_buffer& buf);
    virtual BCP_user_data*
    unpack_user_data(BCP_buffer& buf);

    /** What is the process id of the current process */
    const BCP_proc_id*
    process_id() const;
    /** Send a message to a particular process */
    void
    send_message(const BCP_proc_id* const target, const BCP_buffer& buf);
    /** Broadcast the message to all processes of the given type */
    void
    broadcast_message(const BCP_process_t proc_type, const BCP_buffer& buf);
    /** Process a message that has been sent by another process' user part to
	this process' user part. */
    virtual void
    process_message(BCP_buffer& buf);

    virtual OsiSolverInterface *
    initialize_solver_interface();

    virtual BCP_solution*
    test_feasibility(const BCP_lp_result& lp_result,
		     const BCP_vec<BCP_var*>& vars,
		     const BCP_vec<BCP_cut*>& cuts);
    BCP_solution* test_feasibility_BB(const BCP_vec<BCP_var*>& vars);
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
				BCP_vec<BCP_lp_branching_object*>& cans);

    BCP_branching_decision bbBranch(OsiBranchingInformation& brInfo,
				    BCP_vec<BCP_lp_branching_object*>& cands);
    BCP_branching_decision hybridBranch();

    virtual void
    set_user_data_for_children(BCP_presolved_lp_brobj* best, 
			       const int selected);

    virtual void
    modify_lp_parameters(OsiSolverInterface* lp, bool in_strong_branching);

private:
    /* There's no totalTime_ and nodeTime_. Look at the top of BM.cpp */
    //   double totalTime_;
    //   double nodeTime_;
    int in_strong;

};

//#############################################################################

#include "BCP_USER.hpp"

class BM_init : public USER_initialize {

public:

    virtual BCP_tm_user * tm_init(BCP_tm_prob& p,
				  const int argnum,
				  const char * const * arglist);

    virtual BCP_lp_user * lp_init(BCP_lp_prob& p);
};

//#############################################################################

class BM_solution : public BCP_solution { 
public:
    double _objective;
    BCP_vec<int> _ind;
    BCP_vec<double> _values;

public:
    BM_solution() : _objective(0), _ind(), _values() {}
    virtual ~BM_solution() {}

    inline virtual double objective_value() const { return _objective; }
    inline void setObjective(double obj) { _objective = obj; }

    void add_entry(int i, double value) {
	_ind.push_back(i);
	_values.push_back(value);
    }
};

#endif
