// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
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
#include "BCP_tm_user.hpp"
#include "BCP_lp_user.hpp"

#include "BM_var.hpp"
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
    
class BM_par {
public:
    enum chr_params {
        //
	BranchOnSos,
        CombinedDistanceAndPriority,
        PureBranchAndBound,
	PrintBranchingInfo,
        LowPriorityImportant,
        end_of_chr_params
    };
    enum int_params {
        //
        NumNlpFailureMax,
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

    virtual void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) {
	BM_pack_var(var, buf);
    }
    virtual BCP_var_algo* unpack_var_algo(BCP_buffer& buf) {
	return BM_unpack_var(buf);
    }

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
		BCP_user_data*& user_data,
		BCP_pricing_status& pricing_status);

    /// Print a feasible solution
    virtual void display_feasible_solution(const BCP_solution* sol);

    void readIpopt();

};

//#############################################################################

#include <OsiAuxInfo.hpp>
#include <OsiCuts.hpp>
#include "BCP_lp_user.hpp"
#include "IpCbcOACutGenerator2.hpp"
#include "BonminAmplInterface.hpp"
class BM_lp : public BCP_lp_user
{
    struct BmSosInfo {
	int num;
	int *priority_order;
	int *lengths;
	int *types; // 0: type 1  ---- 1: type 2
	int *priorities;
	int **indices;
	double **weights;
	BmSosInfo() : num(0) {}
	~BmSosInfo();
	void setFrom(const TMINLP::SosInfo * sos);
	void shuffle();
    };
		
    BCP_string ipopt_file_content;
    BCP_string nl_file_content;
    BCP_parameter_set<BM_par> par;

    OsiBabSolver babSolver_;
    BonminAmplInterface nlp;
    BmSosInfo sos;

    double lower_bound_;
    double* primal_solution_;

    /** A counter for how many times in a row did the NLP code fail. When the
	NLP fails we branch; hopefully it'll be OK in the children. If it
	fails too many times in a row then we fathom the node: it's hopelessly
	difficult. */
    int numNlpFailed_;

    /** If pure branch and bound is done then for each fractional variable
	that ought to be integer we multiply it's distance from integrality
	with it's priority and choose the var with the highest value */
    double* branching_priority_;

    IpCbcOACutGenerator2* feasChecker_;
    OsiCuts cuts_;

public:
    BM_lp();
    virtual ~BM_lp();

    inline int& numNlpFailed() {
	return (dynamic_cast<BM_node*>(get_user_data()))->numNlpFailed_;
    }

    virtual void
    unpack_module_data(BCP_buffer& buf);

    virtual void pack_var_algo(const BCP_var_algo* var, BCP_buffer& buf) {
	BM_pack_var(var, buf);
    }
    virtual BCP_var_algo* unpack_var_algo(BCP_buffer& buf) {
	return BM_unpack_var(buf);
    }

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

    virtual OsiSolverInterface *
    initialize_solver_interface();

    virtual BCP_solution*
    test_feasibility(const BCP_lp_result& lp_result,
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
    virtual void
    vars_to_cols(const BCP_vec<BCP_cut*>& cuts,
		 BCP_vec<BCP_var*>& vars,
		 BCP_vec<BCP_col*>& cols,
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
