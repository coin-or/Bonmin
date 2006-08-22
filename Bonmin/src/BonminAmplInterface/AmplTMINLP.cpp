// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 12/01/2004
#include "IpBlas.hpp"

#include "AmplTNLP.hpp"
#include "AmplTMINLP.hpp"
#include <iostream>

#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"
#include "CoinHelperFunctions.hpp"

namespace ampl_utils
{
  void sos_kludge(int nsos, int *sosbeg, double *sosref);
}
namespace Ipopt
{

  AmplTMINLP::AmplTMINLP()
      :
      TMINLP(),
      ampl_tnlp_(NULL),
      branch_(),
      sos_()
  {}


  AmplTMINLP::AmplTMINLP(const SmartPtr<const Journalist>& jnlst,
      const SmartPtr<OptionsList> options,
      char**& argv,
      AmplSuffixHandler* suffix_handler /*=NULL*/,
      const std::string& appName,
      std::string* nl_file_content /* = NULL */)
      :
      TMINLP(),
      ampl_tnlp_(NULL),
      branch_(),
      sos_(),
      suffix_handler_(NULL)
  {
    Initialize(jnlst, options, argv, suffix_handler, appName, nl_file_content);
  }

  void
  AmplTMINLP::Initialize(const SmartPtr<const Journalist>& jnlst,
      const SmartPtr<OptionsList> options,
      char**& argv,
      AmplSuffixHandler* suffix_handler /*=NULL*/,
      const std::string& appName,
      std::string* nl_file_content /* = NULL */)
  {


    if(suffix_handler==NULL)
      suffix_handler_ = suffix_handler = new AmplSuffixHandler();

    // Add the suffix handler for scaling
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
  
    // priority suffix
    suffix_handler->AddAvailableSuffix("priority", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("direction", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("downPseudocost", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("upPseudocost", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);



    // sos suffixes
    suffix_handler->AddAvailableSuffix("ref", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sos",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sos",AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sosno",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sosref",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sstatus", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sstatus", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);
    

    SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();
    fillAmplOptionList(GetRawPtr(ampl_options_list));
    fillApplicationOptions(GetRawPtr(ampl_options_list) );
    std::string options_id = appName + "_options";
    ampl_tnlp_ = new AmplTNLP(jnlst, options, argv, suffix_handler, true,
         ampl_options_list, options_id.c_str(),
        appName.c_str(), appName.c_str(), nl_file_content);
    /* Read suffixes */
    read_priorities();
    read_sos();
  }

  AmplTMINLP::~AmplTMINLP()
  {}

  void
  AmplTMINLP::read_priorities()
  {
    int numcols, m, dummy1, dummy2;
    TNLP::IndexStyleEnum index_style;
    ampl_tnlp_->get_nlp_info(numcols, m, dummy1, dummy2, index_style);

    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);

    const Index* pri = suffix_handler->GetIntegerSuffixValues("priority", AmplSuffixHandler::Variable_Source);
    const Index* brac = suffix_handler->GetIntegerSuffixValues("direction", AmplSuffixHandler::Variable_Source);
    const Number* upPs = suffix_handler->GetNumberSuffixValues("upPseudocost", AmplSuffixHandler::Variable_Source);
    const Number* dwPs = suffix_handler->GetNumberSuffixValues("downPseudocost", AmplSuffixHandler::Variable_Source);


    branch_.gutsOfDestructor();
    branch_.size = numcols;
    if(pri) {
      branch_.priorities = new int[numcols];
      for(int i = 0 ; i < numcols ; i++) {
        branch_.priorities [i] = -pri[i] + 9999;
      }
    }
    if(brac) {
      branch_.branchingDirections = CoinCopyOfArray(brac,numcols);
    }
    if(upPs && !dwPs) dwPs = upPs;
    else if(dwPs && !upPs) upPs = dwPs;
  
    if(upPs) {
      branch_.upPsCosts = CoinCopyOfArray(upPs,numcols);
    }
    if(dwPs) {
      branch_.downPsCosts = CoinCopyOfArray(dwPs,numcols);
    }
  } 

  void
  AmplTMINLP::read_sos()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();

    int i = ASL_suf_sos_explict_free;
    int copri[2], **p_sospri;
    copri[0] = 0;
    copri[1] = 0;
    int * starts = NULL;
    int * indices = NULL;
    char * types = NULL;
    double * weights = NULL;
    int * priorities = NULL;
    p_sospri = &priorities;

    sos_.gutsOfDestructor();

    sos_.num = suf_sos(i, &sos_.numNz, &types, p_sospri, copri,
        &starts, &indices, &weights);
    if (sos_.num) {
      //Copy sos information
      sos_.priorities = CoinCopyOfArray(priorities,sos_.num);
      sos_.starts = CoinCopyOfArray(starts, sos_.num + 1);
      sos_.indices = CoinCopyOfArray(indices, sos_.numNz);
      sos_.types = CoinCopyOfArray(types, sos_.num);
      sos_.weights = CoinCopyOfArray(weights, sos_.numNz);

      ampl_utils::sos_kludge(sos_.num, sos_.starts, sos_.weights);
      for (int ii=0;ii<sos_.num;ii++) {
        int ichar = sos_.types[ii];
        if(ichar != '1') {
          std::cerr<<"Unsuported type of sos constraint: "<<sos_.types[ii]<<std::endl;
          throw;
        }
        sos_.types[ii]= 1;
      }
    }
  }

  bool AmplTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    return ampl_tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  }

  bool AmplTMINLP::get_var_types(Index n, VariableType* var_types)
  {
    // The variables are sorted by type in AMPL, so all we need to
    // know are the counts of each type.


    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);
    int totalNumberOfNonContinuous = 0;
    //begin with variables non-linear both in the objective function and in some constraints
    // The first ones are continuous
    Index start = 0;
    Index end = n_non_linear_b - n_non_linear_bi;
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_bi;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }

    //next variables non-linear in some constraints
    // The first ones are continuous
    start = end;
    end = max(start,end + n_non_linear_c - n_non_linear_ci - n_non_linear_b);
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_ci;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }

    //next variables non-linear in the objective function
    // The first ones are continuous
    start = end;
    end = max(start,end + n_non_linear_o - max(n_non_linear_b, n_non_linear_c) - n_non_linear_oi);
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_oi;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }

    //At last the linear variables
    // The first ones are continuous
    start = end;
    end = n - n_binaries - n_integers;
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are binaries
    start = end;
    end = start + n_binaries;
    for (int i=start; i < end; i++) {
      var_types[i] = BINARY;
      totalNumberOfNonContinuous++;
    }

    // The third ones are integers
    start = end;
    end = start + n_integers;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }
    //    std::cout<<"Number of integer and binaries : "<<totalNumberOfNonContinuous<<std::endl;
    return true;
  }

  bool AmplTMINLP::get_constraints_types(Index n, ConstraintType* const_types)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    //check that n is good
    DBG_ASSERT(n == n_con);
    // check that there are no network constraints
    DBG_ASSERT(nlnc == 0 && lnc == 0);
    //the first nlc constraints are non linear the rest is linear
    int i;
    for(i = 0 ; i < nlc ; i++)
      const_types[i]=NON_LINEAR;
    // the rest is linear
    for(; i < n_con ; i++)
      const_types[i]=LINEAR;
    return true;
  }

  bool AmplTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
      Index m, Number* g_l, Number* g_u)
  {
    return ampl_tnlp_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  }

  bool AmplTMINLP::get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda, Number* lambda)
  {
    return ampl_tnlp_->get_starting_point(n, init_x, x,
        init_z, z_L, z_U,
        m, init_lambda, lambda);
  }

  bool AmplTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    return ampl_tnlp_->eval_f(n, x, new_x, obj_value);
  }

  bool AmplTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    return ampl_tnlp_->eval_grad_f(n, x, new_x, grad_f);
  }

  bool AmplTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    return ampl_tnlp_->eval_g(n, x, new_x, m, g);
  }

  bool AmplTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
      Index m, Index nele_jac, Index* iRow,
      Index *jCol, Number* values)
  {
    return ampl_tnlp_->eval_jac_g(n, x, new_x,
        m, nele_jac, iRow, jCol,
        values);
  }

  bool AmplTMINLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess, Index* iRow,
      Index* jCol, Number* values)
  {
    return ampl_tnlp_->eval_h(n, x, new_x, obj_factor,
        m, lambda, new_lambda, nele_hess, iRow,
        jCol, values);
  }

  void AmplTMINLP::finalize_solution(SolverReturn status,
      Index n, const Number* x, const Number* z_L, const Number* z_U,
      Index m, const Number* g, const Number* lambda,
      Number obj_value)
  {
    // Not sure if ampl require a different form of solution file
    // for MINLPs - we may have to write a different solution file here instead of
    // passing this back to ampl.
    ampl_tnlp_->finalize_solution(status,
        n, x, z_L, z_U,
        m, g, lambda,
        obj_value);

    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    solve_result_num = 0;
  }

  void AmplTMINLP::write_solution(const std::string & message, const Number *x_sol, const Number * lambda_sol)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    ;
    DBG_ASSERT(asl);
    //    DBG_ASSERT(x_sol);

    // We need to copy the message into a non-const char array to make
    // it work with the AMPL C function.
    char* cmessage = new char[message.length()+1];
    strcpy(cmessage, message.c_str());

    // In order to avoid casting into non-const, we copy the data here...
    double* x_sol_copy = NULL;
    if (x_sol) {
      x_sol_copy = new double[n_var];
      for (int i=0; i<n_var; i++) {
        x_sol_copy[i] = x_sol[i];
      }
    }
    double* lambda_sol_copy = NULL;
    if (lambda_sol_copy) {
      lambda_sol_copy = new double[n_con];
      for (int i=0; i<n_con; i++) {
        lambda_sol_copy[i] = lambda_sol[i];
      }
    }
    write_sol(cmessage, x_sol_copy, lambda_sol_copy, NULL);

    delete [] x_sol_copy;
    delete [] lambda_sol_copy;
    delete [] cmessage;

  }


  /** This methods gives the linear part of the objective function */
  void AmplTMINLP::getLinearPartOfObjective(double * obj)
  {
    Index n, m, nnz_jac_g, nnz_h_lag;
    TNLP::IndexStyleEnum index_style = TNLP::FORTRAN_STYLE;
    get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
    eval_grad_f(n, NULL, 0,obj);
    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);

    // Now get the coefficients of variables wich are linear in the objective
    // The first ones are not
    Index start = 0;
    Index end = n_non_linear_b;
    for (int i=start; i<end; i++) {
      obj[i] = 0.;
    }

    //next variables should be linear in the objective just skip them
    // (from current end to (end + n_non_linear_c - n_non_linear_ci - n_non_linear_b;)


    // Those are nonlinear in the objective
    start = end + n_non_linear_c;
    end = start + n_non_linear_o;
    for (int i=start; i < end; i++) {
      obj[i]=0.;
    }
    //The rest is linear keep the values of the gradient
  }


  void
  AmplTMINLP::fillAmplOptionList(AmplOptionsList* amplOptList)
  {
    amplOptList->AddAmplOption("bonmin.algorithm","bonmin.algorithm",
        AmplOptionsList::String_Option,
        "Choice of the algorithm.");

    amplOptList->AddAmplOption("bonmin.bb_log_level","bonmin.bb_log_level",
        AmplOptionsList::Integer_Option,
        "specify BB log level");

    amplOptList->AddAmplOption("bonmin.lp_log_level","bonmin.lp_log_level",
        AmplOptionsList::Integer_Option,
        "specify sub-LP log level");

    amplOptList->AddAmplOption("bonmin.milp_log_level","bonmin.milp_log_level",
        AmplOptionsList::Integer_Option,
        "specify sub-MILP log level");

    amplOptList->AddAmplOption("bonmin.oa_log_level","bonmin.oa_log_level",
        AmplOptionsList::Integer_Option,
        "specify OA log level");

    amplOptList->AddAmplOption("bonmin.oa_log_frequency","bonmin.oa_log_frequency",
        AmplOptionsList::Number_Option,
        "specify OA log frequency");

    amplOptList->AddAmplOption("bonmin.nlp_log_level","bonmin.nlp_log_level",
        AmplOptionsList::Integer_Option,
        "specify sub-NLP log level");

    amplOptList->AddAmplOption("bonmin.print_user_options","bonmin.print_user_options",
        AmplOptionsList::String_Option,
        "print options list");

    amplOptList->AddAmplOption("bonmin.bb_log_interval","bonmin.bb_log_interval",
        AmplOptionsList::Integer_Option,
        "Interval at which bound output is given");

    amplOptList->AddAmplOption("bonmin.allowable_gap","bonmin.allowable_gap",
        AmplOptionsList::Number_Option,
        "Specify allowable absolute gap");

    amplOptList->AddAmplOption("bonmin.allowable_fraction_gap","bonmin.allowable_fraction_gap",
        AmplOptionsList::Number_Option,
        "Specify allowable relative gap");

    amplOptList->AddAmplOption("bonmin.cutoff_decr","bonmin.cutoff_decr",
        AmplOptionsList::Number_Option,
        "Specify cutoff decrement");

    amplOptList->AddAmplOption("bonmin.cutoff","bonmin.cutoff",
        AmplOptionsList::Number_Option,
        "Specify cutoff");

    amplOptList->AddAmplOption("bonmin.nodeselect_stra","bonmin.nodeselect_stra",
        AmplOptionsList::String_Option,
        "Choose the node selection strategy");


    amplOptList->AddAmplOption("bonmin.number_strong_branch", "bonmin.number_strong_branch",
        AmplOptionsList::Integer_Option,
        "Chooes number of variable for strong branching");

    amplOptList->AddAmplOption("bonmin.number_before_trust", "bonmin.number_before_trust",
        AmplOptionsList::Integer_Option,
        "Set number of branches on a variable before its pseudo-costs are to be believed");

    amplOptList->AddAmplOption("bonmin.time_limit", "bonmin.time_limit",
        AmplOptionsList::Number_Option,
        "Set maximum computation time for Algorithm");

    amplOptList->AddAmplOption("bonmin.node_limit","bonmin.node_limit",
        AmplOptionsList::Integer_Option,
        "Set maximum number of nodes explored");

    amplOptList->AddAmplOption("bonmin.integer_tolerance", "bonmin.integer_tolerance",
        AmplOptionsList::Number_Option,
        "Set integer tolerance");

    amplOptList->AddAmplOption("bonmin.warm_start","bonmin.warm_start",
        AmplOptionsList::String_Option,
        "Set warm start method");

    amplOptList->AddAmplOption("bonmin.sos_constraints","bonmin.sos_constraints",
        AmplOptionsList::String_Option,
        "Disable SOS contraints");

    amplOptList->AddAmplOption("bonmin.max_random_point_radius",
        "bonmin.max_random_point_radius",
        AmplOptionsList::Number_Option,
        "Set max value for a random point");

    amplOptList->AddAmplOption("bonmin.max_consecutive_failures",
        "bonmin.max_consecutive_failures",
        AmplOptionsList::Integer_Option,
        "Number of consecutive unsolved problems before aborting.");

    amplOptList->AddAmplOption("bonmin.num_iterations_suspect",
        "bonmin.num_iterations_suspect",
        AmplOptionsList::Integer_Option,
        "Number of iteration over which a node is considered suspect");

    amplOptList->AddAmplOption("bonmin.nlp_failure_behavior",
        "bonmin.nlp_failure_behavior",
        AmplOptionsList::String_Option,
        "Set the behavior when the nlp fails.");

    amplOptList->AddAmplOption("bonmin.num_retry_unsolved_random_point",
        "bonmin.num_retry_unsolved_random_point",
        AmplOptionsList::Integer_Option,
        "Number of tries to resolve a failed NLP with a random starting point");

    amplOptList->AddAmplOption("bonmin.max_consecutive_infeasible",
        "bonmin.max_consecutive_infeasible",
        AmplOptionsList::Integer_Option,
        "Number of consecutive infeasible problems before continuing a"
        " branch.");

    amplOptList->AddAmplOption("bonmin.num_resolve_at_root", "bonmin.num_resolve_at_root",
        AmplOptionsList::Integer_Option,
        "Number of tries to resolve the root node with different starting point (only usefull in non-convex).");

    amplOptList->AddAmplOption("bonmin.num_resolve_at_node","bonmin.num_resolve_at_node",
        AmplOptionsList::Integer_Option,
        "Number of tries to resolve a non root node with different starting point (only usefull in non-convex).");


    amplOptList->AddAmplOption("bonmin.nlp_solve_frequency","bonmin.nlp_solve_frequency",
        AmplOptionsList::Integer_Option,
        "Specify the frequency at which nlp relaxations are solved in hybrid.");

    amplOptList->AddAmplOption("bonmin.oa_dec_time_limit", "bonmin.oa_dec_time_limit",
        AmplOptionsList::Number_Option,
        "Specify the maximum amount of time spent in OA decomposition iteratrions.");

    amplOptList->AddAmplOption("bonmin.tiny_element","bonmin.tiny_element",
        AmplOptionsList::Number_Option,
        "Value for tiny element in OA cut");

    amplOptList->AddAmplOption("bonmin.very_tiny_element","bonmin.very_tiny_element",
        AmplOptionsList::Number_Option,
        "Value for very tiny element in OA cut");

    amplOptList->AddAmplOption("bonmin.milp_subsolver", "bonmin.milp_subsolver",
        AmplOptionsList::String_Option,
        "Choose the subsolver to solve MILPs sub-problems in OA decompositions.");

    amplOptList->AddAmplOption("bonmin.Gomory_cuts", "bonmin.Gomory_cuts",
        AmplOptionsList::Integer_Option,
        "Frequency for Generating Gomory cuts in branch-and-cut.");

    amplOptList->AddAmplOption("bonmin.probing_cuts", "bonmin.probing_cuts",
        AmplOptionsList::Integer_Option,
        "Frequency for Generating probing cuts in branch-and-cut");

    amplOptList->AddAmplOption("bonmin.cover_cuts", "bonmin.cover_cuts",
        AmplOptionsList::Integer_Option,
        "Frequency for Generating cover cuts in branch-and-cut");


    amplOptList->AddAmplOption("bonmin.mir_cuts", "bonmin.mir_cuts",
        AmplOptionsList::Integer_Option,
        "Frequency for Generating MIR cuts in branch-and-cut");

  }
} // namespace Ipopt
