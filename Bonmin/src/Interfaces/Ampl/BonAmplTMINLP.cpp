// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
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
#include "BonAmplTMINLP.hpp"
#include <iostream>

#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

#include <fstream>

#include "CoinHelperFunctions.hpp"
#include "BonExitCodes.hpp"
namespace ampl_utils
{
  void sos_kludge(int nsos, int *sosbeg, double *sosref);
}
namespace Bonmin
{
  
  AmplTMINLP::AmplTMINLP()
  :
  TMINLP(),
  upperBoundingObj_(-1),
  ampl_tnlp_(NULL),
  branch_(),
  sos_(),
  suffix_handler_(NULL),
  constraintsConvexities_(NULL),
  nonConvexConstraintsAndRelaxations_(NULL),
  simpleConcaves_(NULL),
  hasLinearObjective_(false)
  {}
  
  
  AmplTMINLP::AmplTMINLP(const SmartPtr<const Journalist>& jnlst,
                         const SmartPtr<OptionsList> options,
                         char**& argv,
                         AmplSuffixHandler* suffix_handler /*=NULL*/,
                         const std::string& appName,
                         std::string* nl_file_content /* = NULL */
                         )
  :
  TMINLP(),
  upperBoundingObj_(-1),
  ampl_tnlp_(NULL),
  branch_(),
  sos_(),
  suffix_handler_(NULL),
  constraintsConvexities_(NULL),
  nonConvexConstraintsAndRelaxations_(NULL),
  simpleConcaves_(NULL),
  hasLinearObjective_(false)
  {
    Initialize(jnlst, options, argv, suffix_handler, appName, nl_file_content);
  }
  
  void
  AmplTMINLP::Initialize(const SmartPtr<const Journalist>& jnlst,
                         const SmartPtr<OptionsList> options,
                         char**& argv,
                         AmplSuffixHandler* suffix_handler /*=NULL*/,
                         const std::string& appName,
                         std::string* nl_file_content /* = NULL */
                         )
  {
    options->GetEnumValue("file_solution",writeAmplSolFile_,"bonmin.");
    jnlst_ = jnlst;

    if (suffix_handler==NULL)
      suffix_handler_ = suffix_handler = new AmplSuffixHandler();
    
    // Add the suffix handler for scaling
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
    
    // Modified for warm-start from AMPL
    suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

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
    

   // For marking convex/nonconvex constraints
   suffix_handler->AddAvailableSuffix("id",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
   suffix_handler->AddAvailableSuffix("primary_var",AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);
    
    // For objectives
    suffix_handler->AddAvailableSuffix("UBObj", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Index_Type);
    
    
    // Perturbation radius
    suffix_handler->AddAvailableSuffix("perturb_radius",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    
    SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();
    fillAmplOptionList(GetRawPtr(ampl_options_list));
    fillApplicationOptions(GetRawPtr(ampl_options_list) );
    std::string options_id = appName + "_options";
    ampl_tnlp_ = new AmplTNLP(jnlst, options, argv, suffix_handler, true,
                              ampl_options_list, options_id.c_str(),
                              appName.c_str(), appName.c_str(), nl_file_content);
    /* Read suffixes */
    read_obj_suffixes();
    read_priorities();
    read_convexities();
    read_sos();


    /* Determine if objective is linear.*/
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
    if(n_non_linear_b == 0 && n_non_linear_o == 0){
        std::cout<<"Problem has a linear objective"<<std::endl;
        hasLinearObjective_ = true;}
  }
  
  AmplTMINLP::~AmplTMINLP()
  {
    delete [] constraintsConvexities_;
    delete [] simpleConcaves_;
    delete [] nonConvexConstraintsAndRelaxations_;
    delete ampl_tnlp_;
  }
  
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
    if (pri) {
      branch_.priorities = new int[numcols];
      for (int i = 0 ; i < numcols ; i++) {
        branch_.priorities [i] = -pri[i] + 9999;
      }
    }
    if (brac) {
      branch_.branchingDirections = CoinCopyOfArray(brac,numcols);
    }
    if (upPs && !dwPs) dwPs = upPs;
    else if (dwPs && !upPs) upPs = dwPs;
    
    if (upPs) {
      branch_.upPsCosts = CoinCopyOfArray(upPs,numcols);
    }
    if (dwPs) {
      branch_.downPsCosts = CoinCopyOfArray(dwPs,numcols);
    }
    
    const double* perturb_radius =
      suffix_handler->GetNumberSuffixValues("perturb_radius", AmplSuffixHandler::Variable_Source);
    perturb_info_.SetPerturbationArray(numcols, perturb_radius);
  }
  
  void
  AmplTMINLP::read_sos()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    
    int i = 0;//ASL_suf_sos_explict_free;
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
        if (ichar != '1' && ichar != '2') {
          std::cerr<<"Unsuported type of sos constraint: "<<sos_.types[ii]<<std::endl;
          throw;
        }
        sos_.types[ii]= 1;
      }
    }
  }
  
  void
  AmplTMINLP::read_obj_suffixes(){
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    DBG_ASSERT(asl);
    //Get the values
    if(n_obj < 2) return;
    
    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);
    
    const Index* UBObj = suffix_handler->GetIntegerSuffixValues("UBObj", AmplSuffixHandler::Objective_Source);
    if(UBObj){
      //get index of lower bounding objective
      int lowerBoundingObj = (UBObj[0] == 1)? 1:0;
      // Pass information to Ipopt
      ampl_tnlp_->set_active_objective(lowerBoundingObj);
      
      //get the index of upper bounding objective
      for(int i = 0; i < n_obj; i++){
        if(UBObj[i]==1){
          if(upperBoundingObj_ != -1){
            jnlst_->Printf(J_ERROR, J_MAIN,
                                   "Too many objectives for upper-bounding");
          }
          upperBoundingObj_ = i;
        }
      }
    }
    else
    {
      ampl_tnlp_->set_active_objective(0);
    }
  }
 /** To store all data stored in the nonconvex suffixes.*/
 struct NonConvexSuff{
    /** Default constructor.*/
    NonConvexSuff():
    cIdx(-1),relIdx(-1),scXIdx(-1),scYIdx(-1){}
    /** Copy constructor.*/
    NonConvexSuff(const NonConvexSuff& other):
    cIdx(other.cIdx), relIdx(other.relIdx),
    scXIdx(other.scXIdx), scYIdx(other.scYIdx){}
    /** Index of the nonconvex constraint.*/
    int cIdx;
    /** Index of its relaxation.*/
    int relIdx;
    /** Index of variable x in a simple concave constraint of type y >= F(x).*/
    int scXIdx;
    /** Index of variable y in a simple concave constraint of type y >= F(x).*/
    int scYIdx;
    };
 void AmplTMINLP::read_convexities(){
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    DBG_ASSERT(asl);

   const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);
   const Index * id = suffix_handler->GetIntegerSuffixValues("id", AmplSuffixHandler::AmplSuffixHandler::Variable_Source);
   const Index * primary_var = suffix_handler->GetIntegerSuffixValues("primary_var", AmplSuffixHandler::AmplSuffixHandler::Constraint_Source);

   if(primary_var!= NULL)
   {
     if(constraintsConvexities_ != NULL){
       delete [] constraintsConvexities_;}
     constraintsConvexities_ = new TMINLP::Convexity[n_con];
     if(id == NULL){
        std::cout<<"Incorrect suffixes description in ampl model. id's are not declared "<<std::endl;
      exit(ERROR_IN_AMPL_SUFFIXES);
	}
     int numberSimpleConcave = 0; 
     std::map<int, int> id_map;

     for(int i = 0 ; i < n_var ; i++){
        id_map[id[i]] = i;}


     for(int i = 0 ; i < n_con ; i++){
        if(primary_var[i] != 0){
            constraintsConvexities_[i] = TMINLP::SimpleConcave;
            numberSimpleConcave++; 
        }
        else constraintsConvexities_[i] = TMINLP::Convex;
     }
    simpleConcaves_ = new SimpleConcaveConstraint[numberSimpleConcave];
    nonConvexConstraintsAndRelaxations_ = new MarkedNonConvex[numberSimpleConcave];
    numberSimpleConcave = 0;
    int * jCol = new int[n_var];
    for(int i = 0 ; i < n_con ; i++){
        if(primary_var[i] != 0){
           nonConvexConstraintsAndRelaxations_[numberSimpleConcave].cIdx = i;
           nonConvexConstraintsAndRelaxations_[numberSimpleConcave].cRelaxIdx = -1;
           simpleConcaves_[numberSimpleConcave].cIdx = i;
           simpleConcaves_[numberSimpleConcave].yIdx = id_map[primary_var[i]];
           
        //Now get gradient of i to get xIdx.
        int nnz;
        int & yIdx = simpleConcaves_[numberSimpleConcave].yIdx;
	   int & xIdx = simpleConcaves_[numberSimpleConcave].xIdx;
        eval_grad_gi(n_var, NULL, false, i, nnz, jCol, NULL);
        if(nnz != 2){//Error in ampl model
        std::cout<<"Incorrect suffixes description in ampl model. Constraint with id "
                 <<id<<" is simple concave and should have only two nonzero elements"<<std::endl;
         exit(ERROR_IN_AMPL_SUFFIXES);
        }
        if(jCol[0] == yIdx){
          xIdx = jCol[1];
        }
        else{
             if(jCol[1] != yIdx){//Error in ampl model
               std::cout<<"Incorrect suffixes description in ampl model. Constraint with id "
                        <<id<<" : variable marked as y does not appear in the constraint."<<std::endl;
               exit(ERROR_IN_AMPL_SUFFIXES);
             }       
          xIdx = jCol[0];
        }
	   numberSimpleConcave++;
     }
   }
   delete [] jCol;
   numberSimpleConcave_ = numberSimpleConcave;
   numberNonConvex_ = numberSimpleConcave;
   }
   
 }
  
  bool AmplTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    return ampl_tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  }
  
  bool AmplTMINLP::get_variables_types(Index n, VariableType* var_types)
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
  
  
  /** Returns the constraint linearity.
  * array should be alocated with length at least n.*/
  bool 
  AmplTMINLP::get_constraints_linearity(Index m, 
                                        Ipopt::TNLP::LinearityType* const_types)
  {
    return ampl_tnlp_->get_constraints_linearity(m, const_types);
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
  
  bool AmplTMINLP::eval_gi(Index n, const Number* x, bool new_x,
                           Index i, Number& gi)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    
    // ignore new_x for now
    xunknown();
    
    fint nerror = 0;
    gi = conival(i, const_cast<real*>(x), &nerror);
    if (nerror!=0) {
      return false;
    }
    else {
      return true;
    }
  }
  
  bool AmplTMINLP::eval_grad_gi(Index n, const Number* x, bool new_x,
                                Index i, Index& nele_grad_gi, Index* jCol,
                                Number* values)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    
    if (jCol) {
      // Only compute the number of nonzeros and the indices
      DBG_ASSERT(!values);
      nele_grad_gi = 0;
      for (cgrad* cg=Cgrad[i]; cg; cg = cg->next) {
        jCol[nele_grad_gi++] = cg->varno;
      }
      return true;
    }
    DBG_ASSERT(values);
    
    // ignore new_x for now
    xunknown();
    
    asl->i.congrd_mode = 1;
    fint nerror = 0;
    congrd(i, const_cast<real*>(x), values, &nerror);
    if (nerror!=0) {
      return false;
    }
    else {
      return true;
    }
  }
  
  void AmplTMINLP::finalize_solution(TMINLP::SolverReturn status,
                                     Index n, const Number* x, Number obj_value)
  {
    std::string message;
    std::string status_str;
    if(status == TMINLP::SUCCESS) {
      status_str = "\t\"Finished\"";
      message = "\nbonmin: Optimal";
    }
    else if(status == TMINLP::INFEASIBLE) {
      status_str = "\t\"Finished\"";
      message = "\nbonmin: Infeasible problem";
    }
    else if(status == TMINLP::LIMIT_EXCEEDED) {
      status_str = "\t\"Not finished\"";
      message = "\n Optimization interupted on limit.";
    }
    else if(status == TMINLP::MINLP_ERROR) {
      status_str = "\t\"Aborted\"";
      message = "\n Error encountered in optimization.";
    }
    if(writeAmplSolFile_)
    {
      write_solution(message, x);
      std::cout<<"\n "<<status_str<<std::endl;    
    }
   else {
      std::cout<<status_str<<message<<std::endl; 
      std::ofstream of("bonmin.sol");
      for(int i = 0 ; i < n ; i++){
         of<<i<<"\t"<<x[i]<<std::endl;
      }
     of<<"-1\n";
  }
} 

  void AmplTMINLP::write_solution(const std::string & message, const Number *x_sol)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    
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
    write_sol(cmessage, x_sol_copy, NULL, NULL);
    
    delete [] x_sol_copy;
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
    
    amplOptList->AddAmplOption("bonmin.node_comparison","bonmin.node_comparison",
                               AmplOptionsList::String_Option,
                               "Choose the node comparison function");
    
    amplOptList->AddAmplOption("bonmin.tree_search_strategy","bonmin.tree_search_strategy",
                               AmplOptionsList::String_Option,
                               "Choose the node selection strategy");
    
    amplOptList->AddAmplOption("bonmin.varselect_stra","bonmin.varselect_stra",
                               AmplOptionsList::String_Option,
                               "Choose the variable selection strategy");
    
    amplOptList->AddAmplOption("bonmin.varselect_stra","bonmin.varselect_stra",
                               AmplOptionsList::String_Option,
                               "Choose the variable selection strategy");
    
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
    
    amplOptList->AddAmplOption("bonmin.solution_limit","bonmin.solution_limit",
                               AmplOptionsList::Integer_Option,
                               "Set maximum of new best integer before aborting.");
    
    amplOptList->AddAmplOption("bonmin.iteration_limit","bonmin.iteration_limit",
                               AmplOptionsList::Integer_Option,
                               "Set cummulated maximum number of iterations in sub-algorithm used for nodes relaxations");
    
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
    
    amplOptList->AddAmplOption("bonmin.filmint_ecp_cuts","bonmin.filmint_ecp_cuts",
                               AmplOptionsList::Integer_Option,
                               "Specify the frequency (in terms of nodes) at which some a la filmint ecp cuts are generated.");
    
    amplOptList->AddAmplOption("bonmin.couenne_ecp_cuts","bonmin.couenne_ecp_cuts",
                               AmplOptionsList::Integer_Option,
                               "Specify the frequency (in terms of nodes) at which couenne cuts are generated.");
    
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
  
  /** This method to returns the value of an alternative objective function for
  upper bounding (if one has been declared by using the prefix UBObj).*/
  bool 
  AmplTMINLP::eval_upper_bound_f(Index n, const Number* x,
                                  Number& obj_value){
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    //xknown(x);    // This tells ampl to use a new x
    fint nerror = -1;
    obj_value = objval(upperBoundingObj_, const_cast<double *>(x), &nerror);
    if(nerror > 0){
      jnlst_->Printf(J_ERROR, J_MAIN,
                     "Error in evaluating upper bounding objecting");
      throw -1;}
    return nerror;

  }
  /** Return the ampl solver object (ASL*) */
  const ASL_pfgh* 
  AmplTMINLP::AmplSolverObject() const
  {
    return ampl_tnlp_->AmplSolverObject();
  }
  
  
} // namespace Bonmin 
