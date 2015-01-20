// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004, 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 12/01/2004


#include "BonTMINLP2TNLP.hpp"
#include "IpBlas.hpp"
#include "IpAlgTypes.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include <climits>
#include <string>
#include <fstream>
#include <sstream>
#include "Ipopt/BonIpoptInteriorWarmStarter.hpp"
#include "OsiBranchingObject.hpp"

using namespace Ipopt;

extern bool BonminAbortAll;
class OsiObject;
namespace Bonmin
{

  TMINLP2TNLP::TMINLP2TNLP(const SmartPtr<TMINLP> tminlp
#ifdef WARM_STARTER
       ,
      const OptionsList& options
#endif
       )
      :
      var_types_(),
      x_l_(),
      x_u_(),
      orig_x_l_(),
      orig_x_u_(),
      g_l_(),
      g_u_(),
      x_init_(),
      duals_init_(NULL),
      x_init_user_(),
      x_sol_(),
      g_sol_(),
      duals_sol_(),
      tminlp_(tminlp),
      nnz_jac_g_(0),
      nnz_h_lag_(0),
      index_style_(TNLP::FORTRAN_STYLE),
      obj_value_(1e100),
      curr_warm_starter_(),
      nlp_lower_bound_inf_(-DBL_MAX),
      nlp_upper_bound_inf_(DBL_MAX),
      warm_start_entire_iterate_(true),
      need_new_warm_starter_(true)
  {
    // read the nlp size and bounds information from
    // the TMINLP and keep an internal copy. This way the
    // caller can modify the bounds that are sent to Ipopt;
    assert(IsValid(tminlp_));
    Index n,m;
    bool retval =
      tminlp_->get_nlp_info(n, m, nnz_jac_g_, nnz_h_lag_, index_style_);

    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "get_nlp_info of TMINLP returns false.");

    // Allocate space for the variable types vector
    var_types_.resize(n);

    // retrieve the variable types
    tminlp_->get_variables_types(n, var_types_());

    // Allocate space for the internal copy of the variable bounds
    x_l_.resize(n);
    x_u_.resize(n);
    orig_x_l_.resize(n);
    orig_x_u_.resize(n);

    g_l_.resize(m);
    g_u_.resize(m);

    // retrieve the variable bounds
    if(m){
      tminlp_->get_bounds_info(n, x_l_(), x_u_(), m, g_l_(), g_u_());
    }
    else {
      tminlp_->get_bounds_info(n, x_l_(), x_u_(), m, NULL, NULL);
    }
    IpBlasDcopy(n, x_l_(), 1, orig_x_l_(), 1);
    IpBlasDcopy(n, x_u_(), 1, orig_x_u_(), 1);


    // Allocate space for the initial point
    x_init_user_.resize(n);
    tminlp_->get_starting_point(n, true, x_init_user_(), false, NULL, NULL,
        m, false, NULL);

#ifdef WARM_STARTER
    // Get values for parameters
    options.GetNumericValue("nlp_lower_bound_inf", nlp_lower_bound_inf_, "");
    options.GetNumericValue("nlp_upper_bound_inf", nlp_upper_bound_inf_, "");
    options.GetBoolValue("warm_start_entire_iterate",
        warm_start_entire_iterate_, "");
#endif
  }

  TMINLP2TNLP::TMINLP2TNLP(const TMINLP2TNLP& other)
    :
    var_types_(),
    x_l_(),
    x_u_(),
    orig_x_l_(),
    orig_x_u_(),
    g_l_(),
    g_u_(),
    x_init_(),
    duals_init_(NULL),
    x_init_user_(),
    x_sol_(),
    g_sol_(),
    duals_sol_(),
    tminlp_(other.tminlp_),
    nnz_jac_g_(other.nnz_jac_g_),
    nnz_h_lag_(other.nnz_h_lag_),
    index_style_(other.index_style_),
    return_status_(other.return_status_),
    obj_value_(other.obj_value_),
    curr_warm_starter_(other.curr_warm_starter_),
    nlp_lower_bound_inf_(other.nlp_lower_bound_inf_),
    nlp_upper_bound_inf_(other.nlp_upper_bound_inf_),
    warm_start_entire_iterate_(other.warm_start_entire_iterate_),
    need_new_warm_starter_(other.need_new_warm_starter_)
  {
    gutsOfCopy(other);
  }

  /** Overloaded Equals Operator */
  TMINLP2TNLP &
  TMINLP2TNLP::operator=(const TMINLP2TNLP& rhs){
    if(this != &rhs){
      tminlp_ = rhs.tminlp_;
      nnz_jac_g_ = rhs.nnz_jac_g_;
      nnz_h_lag_ = rhs.nnz_h_lag_;
      index_style_ = rhs.index_style_;
      return_status_ = rhs.return_status_;
      obj_value_ = rhs.obj_value_;
      curr_warm_starter_ = rhs.curr_warm_starter_;
      nlp_lower_bound_inf_ = rhs.nlp_lower_bound_inf_;
      nlp_upper_bound_inf_ = rhs.nlp_upper_bound_inf_;
      warm_start_entire_iterate_ = rhs.warm_start_entire_iterate_;
      need_new_warm_starter_ = rhs.need_new_warm_starter_;
  
      gutsOfDelete(); 
      gutsOfCopy(rhs);

    }
    return (*this);
  }

  TMINLP2TNLP::~TMINLP2TNLP()
  {
    gutsOfDelete();
  }

  void
  TMINLP2TNLP::gutsOfDelete(){
  }

  /** Copies all the arrays. 
      \warning this and other should be two instances of the same problem
      \warning AW: I am trying to mimic a copy construction for Cbc
      use with great care not safe.
  */
  void
  TMINLP2TNLP::gutsOfCopy(const TMINLP2TNLP& other)
  {
    Index n = other.num_variables();
    Index m = other.num_constraints();

    if(n > 0){//Copies all the arrays in n_
      var_types_ = other.var_types_;

      x_l_.resize(n);
      x_u_.resize(n); // Those are copied in copyUserModification
      IpBlasDcopy(n, other.x_l_(), 1, x_l_(), 1);
      IpBlasDcopy(n, other.x_u_(), 1, x_u_(), 1);

      orig_x_l_.resize(n);
      orig_x_u_.resize(n);
      IpBlasDcopy(n, other.orig_x_l_(), 1, orig_x_l_(), 1);
      IpBlasDcopy(n, other.orig_x_u_(), 1, orig_x_u_(), 1);
      x_init_user_.resize(n);
      IpBlasDcopy(n, other.x_init_user_(), 1, x_init_user_(), 1);
      if(!other.x_sol_.empty()) {
        Set_x_sol(n,other.x_sol_());
      }
   }

  if(!other.g_l_.empty()){
      const size_t& size = other.g_l_.size();
      g_l_.resize(size);
      g_u_.resize(size);
   }

   if(m > 0){//Copies all the arrays in m_
      IpBlasDcopy(m, other.g_l_(), 1, g_l_(), 1);
      IpBlasDcopy(m, other.g_u_(), 1, g_u_(), 1);
      if(!other.g_sol_.empty()) {
        g_sol_.resize(m);
        IpBlasDcopy(m, other.g_sol_(), 1, g_sol_(), 1);
      }
    }


      x_init_ = other.x_init_;

      if(other.duals_init_) {
        duals_init_ = x_init_() + n;
      }
      else
        duals_init_ = NULL;


    if(!other.duals_sol_.empty()) {
      duals_sol_.resize(m + 2*n);
      IpBlasDcopy((int) duals_sol_.size(), other.duals_sol_(), 1, duals_sol_(), 1);
    }

}

  void TMINLP2TNLP::SetVariablesBounds(Index n,
                                       const Number * x_l,
                                       const Number * x_u)
  {
    assert(n==num_variables());
    IpBlasDcopy(n, x_l, 1, x_l_(), 1);
    IpBlasDcopy(n, x_u, 1, x_u_(), 1);
  }

   void TMINLP2TNLP::SetVariablesLowerBounds(Index n,
                                       const Number * x_l)
  {
    assert(n==num_variables());
    IpBlasDcopy(n, x_l, 1, x_l_(), 1);
  }

   void TMINLP2TNLP::SetVariablesUpperBounds(Index n,
                                       const Number * x_u)
  {
    assert(n==num_variables());
    IpBlasDcopy(n, x_u, 1, x_u_(), 1);
  }

  void TMINLP2TNLP::SetVariableBounds(Index var_no, Number x_l, Number x_u)
  {
    assert(var_no >= 0 && var_no < num_variables());
    x_l_[var_no] = x_l;
    x_u_[var_no] = x_u;
  }

  void TMINLP2TNLP::SetVariableLowerBound(Index var_no, Number x_l)
  {
    assert(var_no >= 0 && var_no < num_variables());
    x_l_[var_no] = x_l;
  }

  void TMINLP2TNLP::SetVariableUpperBound(Index var_no, Number x_u)
  {
    assert(var_no >= 0 && var_no < num_variables());
    x_u_[var_no] = x_u;
  }

  void TMINLP2TNLP::resetStartingPoint()
  {
    curr_warm_starter_ = NULL;
    x_init_.clear();
  }

  void TMINLP2TNLP::setxInit(Index n,const Number* x_init)
  {
    assert(n == num_variables());
    if((int)x_init_.size() < n)
      x_init_.resize(n);
    IpBlasDcopy(n, x_init, 1, x_init_(), 1);
  }

  void TMINLP2TNLP::setDualsInit(Index m, const Number* duals_init)
  {
    assert(m == num_variables() * 2 + num_constraints() );
    x_init_.resize(num_variables() * 3 + num_constraints(), 0.);
    duals_init_ = x_init_() + num_variables();

    if(m >0)
      IpBlasDcopy(m, duals_init, 1, duals_init_, 1);

  }

  /** Set the contiuous solution */
  void TMINLP2TNLP::Set_x_sol(Index n, const Number* x_sol)
  {
    assert(n == num_variables());
    if (x_sol_.empty()) {
      x_sol_.resize(n);
    }
    assert(n == (int) x_sol_.size());
    IpBlasDcopy(n, x_sol, 1, x_sol_(), 1);
  }

  /** Set the contiuous dual solution */
  void TMINLP2TNLP::Set_dual_sol(Index n, const Number* dual_sol)
  {
    assert(n == num_variables() *2 + num_constraints());
    if (duals_sol_.empty()) {
      duals_sol_.resize(n);
    }
    assert(n == (int) duals_sol_.size());
    IpBlasDcopy(n, dual_sol, 1, duals_sol_(), 1);
  }

  /** Change the type of the variable */
  void TMINLP2TNLP::SetVariableType(Index n, TMINLP::VariableType type)
  {
    assert(n >= 0 && n < num_variables());
    var_types_[n] = type;
  }

  bool TMINLP2TNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
      Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    n = num_variables();
    m = num_constraints();
    nnz_jac_g = nnz_jac_g_;
    nnz_h_lag = nnz_h_lag_;
    index_style = index_style_;
     //printf("Been there and said %i\n", nnz_jac_g_);
    return true;
  }

  bool TMINLP2TNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
      Index m, Number* g_l, Number* g_u)
  {
    assert(n==num_variables());
    assert(m==num_constraints());
    IpBlasDcopy(n, x_l_(), 1, x_l, 1);
    IpBlasDcopy(n, x_u_(), 1, x_u, 1);
    if (m > 0){
      IpBlasDcopy(m, g_l_(), 1, g_l, 1);
      IpBlasDcopy(m, g_u_(), 1, g_u, 1);
    }
    return true;
  }

  bool TMINLP2TNLP::get_starting_point(Index n, bool init_x, Number* x,
      bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda,
      Number* lambda)
  {
    assert(m==num_constraints());
    assert(n==num_variables());
#if 0
    x_init_.resize(3*n + m, 0.);
    duals_init_ = x_init_() + n;
#endif
    if (init_x == true) {
      if(x_init_.empty()){
        assert(x_init_user_.size() >= n);
        IpBlasDcopy(n, x_init_user_(), 1, x, 1);
      }
      else
        IpBlasDcopy(n, x_init_(), 1, x, 1);
    }
    if (init_z == true) {
      if(duals_init_ == NULL)
        return false;
      assert(x_init_.size() == 3*n + m && duals_init_ == x_init_() + n); 
      IpBlasDcopy(n, duals_init_, 1, z_L, 1);
      IpBlasDcopy(n, duals_init_ + n, 1, z_U, 1);

    }
    if(init_lambda == true) {
      if(duals_init_ == NULL)
        return false;
      assert(x_init_.size() == 3*n + m && duals_init_ == x_init_() + n); 
      if(m > 0)
        IpBlasDcopy(m, duals_init_ + 2*n , 1, lambda, 1);
    }

    need_new_warm_starter_ = true;
    return true;
  }

  bool TMINLP2TNLP::get_warm_start_iterate(IteratesVector& warm_start_iterate)
  {
    if (IsNull(curr_warm_starter_)) {
      return false;
    }

    bool retval = curr_warm_starter_->WarmStartIterate(num_variables(), x_l_(), x_u_(),
        warm_start_iterate);

    need_new_warm_starter_ = true;
    return retval;
  }

  bool TMINLP2TNLP::eval_f(Index n, const Number* x, bool new_x,
      Number& obj_value)
  {
    return tminlp_->eval_f(n, x, new_x, obj_value);
  }

  bool TMINLP2TNLP::eval_grad_f(Index n, const Number* x, bool new_x,
      Number* grad_f)
  {
    grad_f[n-1] = 0;
    return tminlp_->eval_grad_f(n, x, new_x, grad_f);
  }

  bool TMINLP2TNLP::eval_g(Index n, const Number* x, bool new_x,
      Index m, Number* g)
  {
    int return_code = tminlp_->eval_g(n, x, new_x, m, g);
    return return_code;
  }

  bool TMINLP2TNLP::eval_jac_g(Index n, const Number* x, bool new_x,
      Index m, Index nele_jac, Index* iRow,
      Index *jCol, Number* values)
  {
    bool return_code =
      tminlp_->eval_jac_g(n, x, new_x, m, nele_jac,
			  iRow, jCol, values);
    if(iRow != NULL){
      Index buf;
      for(Index k = 0; k < nele_jac ; k++){
        buf = iRow[k];
        iRow[k] = -1;
        iRow[k] = buf;
      }
    }
    return return_code;
  }

  bool TMINLP2TNLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess,
      Index* iRow, Index* jCol, Number* values)
  {
    return tminlp_->eval_h(n, x, new_x, obj_factor, m, lambda,
        new_lambda, nele_hess,
        iRow, jCol, values);
  }


  bool TMINLP2TNLP::eval_gi(Index n, const Number* x, bool new_x,
                           Index i, Number& gi)
  {
    return tminlp_->eval_gi(n, x, new_x, i, gi);
  }
  
  bool TMINLP2TNLP::eval_grad_gi(Index n, const Number* x, bool new_x,
                                Index i, Index& nele_grad_gi, Index* jCol,
                                Number* values)
  {
    return tminlp_->eval_grad_gi(n, x, new_x, i, nele_grad_gi, jCol, values);
  }

  void TMINLP2TNLP::finalize_solution(SolverReturn status,
      Index n, const Number* x, const Number* z_L, const Number* z_U,
      Index m, const Number* g, const Number* lambda,
      Number obj_value,
      const IpoptData* ip_data,
      IpoptCalculatedQuantities* ip_cq)
  {
    assert(n == (Index) num_variables());
    assert(m == (Index) num_constraints());
    x_sol_.resize(n);
    
    IpBlasDcopy(n, x, 1, x_sol_(), 1);
    
    if(m > 0){
    g_sol_.resize(m);
    IpBlasDcopy(m, g, 1, g_sol_(), 1);
    }
    duals_sol_.resize(m + 2*n);
    if(lambda){
      if(m > 0)
      IpBlasDcopy(m, lambda, 1, duals_sol_() + 2*n, 1);
      
      IpBlasDcopy(n, z_L, 1 , duals_sol_() , 1);
      IpBlasDcopy(n, z_U, 1 , duals_sol_() + n, 1);
    }

    return_status_ = status;
    obj_value_ = obj_value;

    if(status == Ipopt::LOCAL_INFEASIBILITY  && ip_cq != NULL){
      obj_value_ = ip_cq->curr_nlp_constraint_violation(NORM_MAX);
    }
    if (IsValid(curr_warm_starter_)) {
      curr_warm_starter_->Finalize();
    }
  }


  bool TMINLP2TNLP::intermediate_callback(AlgorithmMode mode,
      Index iter, Number obj_value,
      Number inf_pr, Number inf_du,
      Number mu, Number d_norm,
      Number regularization_size,
      Number alpha_du, Number alpha_pr,
      Index ls_trials,
      const IpoptData* ip_data,
      IpoptCalculatedQuantities* ip_cq)
  {
    if (BonminAbortAll) return false;
#if WARM_STARTER
    // If we don't have this swtiched on, we assume that also the
    // "warm_start" option for bonmin is set not to refer to the
    // interior warm start object
    if (!warm_start_entire_iterate_) {
      return true;
    }
    if (need_new_warm_starter_) {
      // Create a new object for later warm start information
      curr_warm_starter_ = new IpoptInteriorWarmStarter(n_, x_l_, x_u_,
          nlp_lower_bound_inf_,
          nlp_upper_bound_inf_,
          warm_start_entire_iterate_);
      need_new_warm_starter_ = false;
    }

    return curr_warm_starter_->UpdateStoredIterates(mode, *ip_data, *ip_cq);
#else 
    return true;
#endif
  }


  /** Procedure to ouptut relevant informations to reproduce a sub-problem.
  Compare the current problem to the problem to solve
  and writes files with bounds which have changed and current starting point.

  */
  void
  TMINLP2TNLP::outputDiffs(const std::string& probName, const std::string * varNames)
  {
    const int &numcols = num_variables();
    const int &numrows = num_constraints();

    const double * currentLower = x_l();
    const double * currentUpper = x_u();

    const double * originalLower = orig_x_l();
    const double * originalUpper = orig_x_u();
    CoinRelFltEq eq;
    std::string fBoundsName = probName;
    std::ostringstream os;
    fBoundsName+=".bounds";
    std::string fModName = probName;
    fModName+= ".mod";
    std::ofstream fBounds;
    std::ofstream fMod;
    bool hasVarNames = 0;

    if(varNames!=NULL )
      hasVarNames=1;
    if(hasVarNames)
      fMod.open(fModName.c_str());
    fBounds.open(fBoundsName.c_str());

    for(int i = 0 ; i < numcols ; i++) {
      if(!eq(currentLower[i],originalLower[i])) {
        if(hasVarNames)
          fMod<<"bounds"<<i<<": "
          <<varNames[i]<<" >= "
          <<currentLower[i]<<";\n";


        fBounds<<"LO"<<"\t"<<i<<"\t"<<currentLower[i]<<std::endl;
      }
      if(!eq(currentUpper[i],originalUpper[i])) {
        if(hasVarNames)
          fMod<<"bounds"<<i<<": "
          <<varNames[i]<<" <= "
          <<currentUpper[i]<<";\n";

        fBounds<<"UP"<<"\t"<<i<<"\t"<<currentUpper[i]<<std::endl;
      }
    }

    //write a file with starting point
    std::string fStartPointName=probName;
    fStartPointName+=".start";

    std::ofstream fStartPoint(fStartPointName.c_str());
    const double * primals = x_init();
    const double * duals = duals_init();
    fStartPoint.precision(17);
    fStartPoint<<numcols<<"\t"<<2*numcols+numrows<<std::endl;
    for(int i = 0 ; i < numcols ; i++)
      fStartPoint<<primals[i]<<std::endl;
    int end = 2*numcols + numrows;
    if(duals) {
      for(int i = 0 ; i < end; i++)
        fStartPoint<<duals[i]<<std::endl;
    }

  }

  /** force solution to be fractionnal.*/
  void
  TMINLP2TNLP::force_fractionnal_sol()
  {
    for(int i=0 ; i < num_variables() ; i++) {
      if( ( var_types_[i] == TMINLP::INTEGER ||
          var_types_[i] == TMINLP::BINARY )&&
          x_l_[i] < x_u_[i] + 0.5)//not fixed
      {
        x_sol_[i] = ceil(x_l_[i]) + 0.5;//make it integer infeasible
      }
    }
  }

  bool 
  TMINLP2TNLP::get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling)
  {
    return tminlp_->get_scaling_parameters(obj_scaling, use_x_scaling, n,
				  x_scaling,
				  use_g_scaling, m, g_scaling);
  }
				  

    /** Method called to check wether a problem has still some variable not fixed. If there are no more
        unfixed vars, checks wether the solution given by the bounds is feasible.*/

    /** @name Methods for setting and getting the warm starter */
  void 
  TMINLP2TNLP::SetWarmStarter(SmartPtr<IpoptInteriorWarmStarter> warm_starter)
    {
      curr_warm_starter_ = warm_starter;
    }
  SmartPtr<IpoptInteriorWarmStarter> 
  TMINLP2TNLP::GetWarmStarter()
    {
      return curr_warm_starter_;
    }


  /** Evaluate the upper bounding function at given point and store the result.*/
  double 
  TMINLP2TNLP::evaluateUpperBoundingFunction(const double * x){
    Number help;
    tminlp_->eval_upper_bound_f(num_variables(), x, help);
    return help;
  }

  double
  TMINLP2TNLP::check_solution(OsiObject ** objects, int nObjects){
    assert(x_sol_.size() == num_variables());
    assert(g_sol_.size() == num_constraints());
    if (objects) {
      for (int i = 0 ; i < nObjects ; i++) {
        OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *>(objects[i]);
        if(obj){
          int colNumber = obj->columnNumber();
          x_sol_[colNumber] = floor(x_sol_[colNumber]+0.5);
        }
      }
    }
    else {
      for (unsigned int i = 0; i < x_sol_.size() ; i++) {
        if (var_types_[i] == TMINLP::INTEGER || var_types_[i] == TMINLP::BINARY) {
          x_sol_[i] = floor(x_sol_[i]+0.5);
        }
      }
    }
    eval_g((int)x_sol_.size(), x_sol_(), true, (int)g_sol_.size(), g_sol_());
    eval_f((int)x_sol_.size(), x_sol_(), false, obj_value_);
    double error = 0;
    for(unsigned int i = 0 ; i < g_sol_.size() ; i++){
      error = std::max(error, std::max(0., g_l_[i] - g_sol_[i]));
      error = std::max(error, std::max(0., - g_u_[i] + g_sol_[i]));
    }
    return error;
  }

}// namespace Bonmin

