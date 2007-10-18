// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
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
#include <climits>
#include <string>
#include <fstream>
#include <sstream>
#include "Ipopt/BonIpoptInteriorWarmStarter.hpp"

extern bool BonminAbortAll;

namespace Bonmin
{

  TMINLP2TNLP::TMINLP2TNLP(const SmartPtr<TMINLP> tminlp
#ifdef WARM_STARTER
       ,
      const OptionsList& options
#endif
       )
      :
      tminlp_(tminlp),
      n_(0),
      m_(0),
      nnz_jac_g_(0),
      nnz_h_lag_(0),
      index_style_(TNLP::FORTRAN_STYLE),
      var_types_(NULL),
      x_l_(NULL),
      x_u_(NULL),
      orig_x_l_(NULL),
      orig_x_u_(NULL),
      g_l_(NULL),
      g_u_(NULL),
      x_init_(NULL),
      duals_init_(NULL),
      capacity_x_init_(0),
      x_init_user_(NULL),
      x_sol_(NULL),
      g_sol_(NULL),
      duals_sol_(NULL),
      obj_value_(1e100),
      curr_warm_starter_(NULL),
      nlp_lower_bound_inf_(-DBL_MAX),
      nlp_upper_bound_inf_(DBL_MAX),
      warm_start_entire_iterate_(true),
      need_new_warm_starter_(true),
      cutsjCol_(NULL),
      cutsiRow_(NULL),
      cutsElems_(NULL),
      cutsLower_(NULL),
      cutsUpper_(NULL),
      nLinearCuts_(0),
      linearCutsNnz_(0),
      linearCutsCapacity_(0),
      linearCutsNnzCapacity_(0) 
  {
    // read the nlp size and bounds information from
    // the TMINLP and keep an internal copy. This way the
    // caller can modify the bounds that are sent to Ipopt;
    assert(IsValid(tminlp_));

    bool retval =
      tminlp_->get_nlp_info(n_, m_, nnz_jac_g_, nnz_h_lag_, index_style_);

    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "get_nlp_info of TMINLP returns false.");

    // Allocate space for the variable types vector
    var_types_ = new TMINLP::VariableType[n_];

    // retrieve the variable types
    tminlp_->get_variables_types(n_, var_types_);

    // Allocate space for the internal copy of the variable bounds
    x_l_ = new Number[n_];
    x_u_ = new Number[n_];
    orig_x_l_ = new Number[n_];
    orig_x_u_ = new Number[n_];

    g_l_ = new Number[m_];
    g_u_ = new Number[m_];

    // retrieve the variable bounds
    tminlp_->get_bounds_info(n_, x_l_, x_u_, m_, g_l_, g_u_);
    IpBlasDcopy(n_, x_l_, 1, orig_x_l_, 1);
    IpBlasDcopy(n_, x_u_, 1, orig_x_u_, 1);

    //    // Check that the bounds make sense compared with the variable type
    //    for (int i=0; i<n_; i++) {
    //      throw_exception_on_bad_variable_bound(i);
    //    }

    // Allocate space for the initial point
    capacity_x_init_ = 3*n_ + 2* m_;//Leave space for adding constraint
    x_init_ = new Number[capacity_x_init_];
    tminlp_->get_starting_point(n_, true, x_init_, false, NULL, NULL,
        m_, false, NULL);
    CoinZeroN(x_init_ + n_ , 2*n_ + 2*m_);
    x_init_user_ = new Number[n_];
    IpBlasDcopy(n_, x_init_, 1, x_init_user_, 1);
    duals_sol_ = NULL;
    duals_init_ = NULL;

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
    tminlp_(other.tminlp_),
    n_(other.n_),
    m_(other.m_),
    nnz_jac_g_(other.nnz_jac_g_),
    nnz_h_lag_(other.nnz_h_lag_),
    index_style_(other.index_style_),
    var_types_(NULL),
    x_l_(NULL),
    x_u_(NULL),
    orig_x_l_(NULL),
    orig_x_u_(NULL),
    g_l_(NULL),
    g_u_(NULL),
    x_init_(NULL),
    duals_init_(NULL),
    capacity_x_init_(other.capacity_x_init_),
    x_init_user_(NULL),
    x_sol_(NULL),
    g_sol_(NULL),
    duals_sol_(NULL),
    return_status_(other.return_status_),
    obj_value_(other.obj_value_),
    curr_warm_starter_(other.curr_warm_starter_),
    nlp_lower_bound_inf_(other.nlp_lower_bound_inf_),
    nlp_upper_bound_inf_(other.nlp_upper_bound_inf_),
    warm_start_entire_iterate_(other.warm_start_entire_iterate_),
    need_new_warm_starter_(other.need_new_warm_starter_),
    cutsjCol_(NULL),
    cutsiRow_(NULL),
    cutsElems_(NULL),
    cutsLower_(NULL),
    cutsUpper_(NULL),
    nLinearCuts_(other.nLinearCuts_),
    linearCutsNnz_(other.linearCutsNnz_),
    linearCutsCapacity_(other.linearCutsCapacity_),
    linearCutsNnzCapacity_(other.linearCutsNnzCapacity_) 
  {
    gutsOfCopy(other);
  }

  /** Overloaded Equals Operator */
  TMINLP2TNLP &
  TMINLP2TNLP::operator=(const TMINLP2TNLP& rhs){
    if(this != &rhs){
      tminlp_ = rhs.tminlp_;
      n_ = rhs.n_;
      m_ = rhs.m_;
      nnz_jac_g_ = rhs.nnz_jac_g_;
      nnz_h_lag_ = rhs.nnz_h_lag_;
      index_style_ = rhs.index_style_;
      capacity_x_init_ = rhs.capacity_x_init_;
      return_status_ = rhs.return_status_;
      obj_value_ = rhs.obj_value_;
      curr_warm_starter_ = rhs.curr_warm_starter_;
      nlp_lower_bound_inf_ = rhs.nlp_lower_bound_inf_;
      nlp_upper_bound_inf_ = rhs.nlp_upper_bound_inf_;
      warm_start_entire_iterate_ = rhs.warm_start_entire_iterate_;
      need_new_warm_starter_ = rhs.need_new_warm_starter_;
      nLinearCuts_ = rhs.nLinearCuts_;
      linearCutsNnz_ = rhs.linearCutsNnz_;
      linearCutsCapacity_ = rhs.linearCutsCapacity_;
      linearCutsNnzCapacity_ = rhs.linearCutsNnzCapacity_;
  
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
    delete [] var_types_;
    var_types_ = NULL;
    delete [] x_l_;
    x_l_ = NULL;
    delete [] x_u_;
    x_u_ = NULL;
    delete [] orig_x_l_;
    orig_x_l_ = NULL;
    delete [] orig_x_u_;
    orig_x_u_ = NULL;
    delete [] g_l_;
    g_l_ = NULL;
    delete [] g_u_;
    g_u_ = NULL;
    delete [] x_init_;
    x_init_ = NULL;
    delete [] x_init_user_;
    x_init_user_ = NULL;
    delete [] x_sol_;
    x_sol_ = NULL;
    delete [] g_sol_;
    g_sol_ = NULL;
    delete [] duals_sol_;
    duals_sol_ = NULL;
    delete [] cutsElems_;
    cutsElems_ = NULL;
    delete [] cutsiRow_;
    cutsiRow_ = NULL;
    delete [] cutsjCol_;
    cutsjCol_ = NULL;
    cutsLower_ = NULL;
    cutsUpper_ = NULL;
  }

  /** Copies all the arrays. 
      \warning this and other should be two instances of the same problem
      \warning AW: I am trying to mimic a copy construction for Cbc
      use with great care not safe.
  */
  void
  TMINLP2TNLP::gutsOfCopy(const TMINLP2TNLP& other)
  {
    assert(n_ == other.n_);
    assert(m_ == other.m_);

    if(n_ > 0){//Copies all the arrays in n_
      var_types_ = new TMINLP::VariableType[n_];
      for (Index i=0; i<n_; i++) {
        var_types_[i] = other.var_types_[i];
      }

      x_l_ = new Number[n_];
      x_u_ = new Number[n_]; // Those are copied in copyUserModification
      IpBlasDcopy(n_, other.x_l_, 1, x_l_, 1);
      IpBlasDcopy(n_, other.x_u_, 1, x_u_, 1);

      orig_x_l_ = new Number[n_];
      orig_x_u_ = new Number[n_];
      IpBlasDcopy(n_, other.orig_x_l_, 1, orig_x_l_, 1);
      IpBlasDcopy(n_, other.orig_x_u_, 1, orig_x_u_, 1);
      x_init_user_ = new Number[n_];
      IpBlasDcopy(n_, other.x_init_user_, 1, x_init_user_, 1);
      if(other.x_sol_ !=NULL) {
        Set_x_sol(n_,other.x_sol_);
      }
   }

  if(m_ + linearCutsCapacity_ > 0){
      int size = m_ + linearCutsCapacity_;
      g_l_ = new Number[size];
      g_u_ = new Number[size];
      cutsLower_ = g_l_ + m_;
      cutsUpper_ = g_u_ + m_;
   }

   if(m_ > 0){//Copies all the arrays in m_
      IpBlasDcopy(m_, other.g_l_, 1, g_l_, 1);
      IpBlasDcopy(m_, other.g_u_, 1, g_u_, 1);
      if(other.g_sol_!=NULL) {
        g_sol_ = new Number [m_];
        IpBlasDcopy(m_, other.g_sol_, 1, g_sol_, 1);
      }
    }


    if(capacity_x_init_ > 0){
      x_init_ = new Number[capacity_x_init_];

      if(other.duals_init_) {
        duals_init_ = &x_init_[n_];
        IpBlasDcopy(capacity_x_init_, other.x_init_, 1, x_init_, 1);
      }
      else
        IpBlasDcopy(n_, other.x_init_, 1, x_init_, 1);
    }


   if(m_ > 0 || n_ > 0){
    if(other.duals_sol_!=NULL) {
      duals_sol_ = new Number[m_ + 2*n_];
      IpBlasDcopy(2*n_+ m_, other.duals_sol_, 1, duals_sol_, 1);
    }
   }

   if(linearCutsNnzCapacity_){
    cutsjCol_ = new int[linearCutsNnzCapacity_];
    cutsiRow_ = new int[linearCutsNnzCapacity_];
    cutsElems_ = new double[linearCutsNnzCapacity_];
    if(linearCutsNnz_ > 0){
      CoinCopyN(other.cutsjCol_, linearCutsNnz_, cutsjCol_);
      CoinCopyN(other.cutsiRow_, linearCutsNnz_, cutsiRow_);
      IpBlasDcopy(linearCutsNnz_, other.cutsElems_, 1, cutsElems_, 1);
    }
  }

  if(nLinearCuts_ > 0){
      IpBlasDcopy(nLinearCuts_, other.cutsLower_, 1, cutsLower_, 1);
      IpBlasDcopy(nLinearCuts_, other.cutsUpper_, 1, cutsUpper_, 1);
  }
}

  void TMINLP2TNLP::SetVariablesBounds(Index n,
                                       const Number * x_l,
                                       const Number * x_u)
  {
    assert(n==n_);
    IpBlasDcopy(n_, x_l, 1, x_l_, 1);
    IpBlasDcopy(n_, x_u, 1, x_u_, 1);
  }

   void TMINLP2TNLP::SetVariablesLowerBounds(Index n,
                                       const Number * x_l)
  {
    assert(n==n_);
    IpBlasDcopy(n_, x_l, 1, x_l_, 1);
  }

   void TMINLP2TNLP::SetVariablesUpperBounds(Index n,
                                       const Number * x_u)
  {
    assert(n==n_);
    IpBlasDcopy(n_, x_u, 1, x_u_, 1);
  }

  void TMINLP2TNLP::SetVariableBounds(Index var_no, Number x_l, Number x_u)
  {
    assert(var_no >= 0 && var_no < n_);
    x_l_[var_no] = x_l;
    x_u_[var_no] = x_u;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetVariableLowerBound(Index var_no, Number x_l)
  {
    assert(var_no >= 0 && var_no < n_);
    x_l_[var_no] = x_l;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetVariableUpperBound(Index var_no, Number x_u)
  {
    assert(var_no >= 0 && var_no < n_);
    x_u_[var_no] = x_u;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetStartingPoint(Index n,const Number* x_init)
  {
    assert(n == n_);
    IpBlasDcopy(n_, x_init, 1, x_init_, 1);
  }

  void TMINLP2TNLP::resetStartingPoint()
  {
    curr_warm_starter_ = NULL;
    IpBlasDcopy(n_, x_init_user_, 1, x_init_, 1);
  }

  void TMINLP2TNLP::setxInit(Index ind, const Number val)
  {
    x_init_[ind] = val;
  }

  void TMINLP2TNLP::setxInit(Index n,const Number* x_init)
  {
    assert(n == n_);
    IpBlasDcopy(n_, x_init, 1, x_init_, 1);
  }

  void TMINLP2TNLP::setDualInit(Index ind, const Number val)
  {
    if(!duals_init_)
      duals_init_ = &x_init_[n_];
    duals_init_[ind] = val;
  }

  void TMINLP2TNLP::setDualsInit(Index m, const Number* duals_init)
  {
    assert(m == m_ + 2*n_ );

    if(!duals_init_)
      duals_init_ = x_init_ + n_;

    IpBlasDcopy(m, duals_init, 1, duals_init_, 1);

  }

  /** Set the contiuous solution */
  void TMINLP2TNLP::Set_x_sol(Index n, const Number* x_sol)
  {
    assert(n == n_);
    if (!x_sol_) {
      x_sol_ = new Number[n];
    }
    IpBlasDcopy(n_, x_sol, 1, x_sol_, 1);
  }

  /** Change the type of the variable */
  void TMINLP2TNLP::SetVariableType(Index n, TMINLP::VariableType type)
  {
    assert(n >= 0 && n < n_);
    var_types_[n] = type;
  }

  bool TMINLP2TNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
      Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    n = n_;
    m = m_ + nLinearCuts_;
    nnz_jac_g = nnz_jac_g_ + linearCutsNnz_;
    nnz_h_lag = nnz_h_lag_;
    index_style = index_style_;
    return true;
  }

  bool TMINLP2TNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
      Index m, Number* g_l, Number* g_u)
  {
    assert(n==n_);
    assert(m==(m_ + (int) nLinearCuts_));
    IpBlasDcopy(n_, x_l_, 1, x_l, 1);
    IpBlasDcopy(n_, x_u_, 1, x_u, 1);
    IpBlasDcopy(m_, g_l_, 1, g_l, 1);
    IpBlasDcopy(m_, g_u_, 1, g_u, 1);
    IpBlasDcopy(nLinearCuts_, cutsLower_, 1, &g_l[m_],1); 
    IpBlasDcopy(nLinearCuts_, cutsUpper_, 1, &g_u[m_],1); 
    return true;
  }

  bool TMINLP2TNLP::get_starting_point(Index n, bool init_x, Number* x,
      bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda,
      Number* lambda)
  {
    assert((int) nLinearCuts_ >= 0);
    assert(m==m_ + (int) nLinearCuts_);
    assert(n==n_);
    if (init_x == true) {
      if(x_init_==NULL)
        return false;
      IpBlasDcopy(n, x_init_, 1, x, 1);
    }
    if (init_z == true) {
      if(duals_init_ == NULL)
        return false;
      IpBlasDcopy(n, duals_init_, 1, z_L, 1);
      IpBlasDcopy(n, duals_init_ + n, 1, z_U, 1);

    }
    if(init_lambda == true) {
      if(duals_init_ == NULL)
        return false;
      int size= 3*n + m_ + nLinearCuts_;
      if(size > capacity_x_init_)
        resizeStartingPoint(0);
      IpBlasDcopy(m_ + nLinearCuts_ , duals_init_ + 2*n , 1, lambda, 1);
    }

    need_new_warm_starter_ = true;
    return true;
  }

  bool TMINLP2TNLP::get_warm_start_iterate(IteratesVector& warm_start_iterate)
  {
    if (IsNull(curr_warm_starter_)) {
      return false;
    }

    bool retval = curr_warm_starter_->WarmStartIterate(n_, x_l_, x_u_,
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

  void TMINLP2TNLP::eval_g_add_linear_cuts(Index n, const Number *x, Index nCuts, Number* g)
  {
    assert(nCuts == nCuts);
    assert(n == n_);
    // Add the linear cuts
    CoinZeroN(g , nLinearCuts_);
    for(unsigned int nnz = 0 ; nnz < linearCutsNnz_ ; nnz++){
      g[cutsiRow_[nnz]] += cutsElems_[nnz] * x[cutsjCol_[nnz]];
    }

  }

  bool TMINLP2TNLP::eval_g(Index n, const Number* x, bool new_x,
      Index m, Number* g)
  {
    assert(m = m_ + nLinearCuts_);
    int return_code = tminlp_->eval_g(n, x, new_x, m_, g);
    if (return_code) {
      eval_g_add_linear_cuts(n, x, nLinearCuts_, 
                             g + m_);
    }
    return return_code;
  }

  void TMINLP2TNLP::eval_jac_g_add_linear_cuts(Index nele_cuts_jac, Index* iRow,
					       Index* jCol, Number* values)
  {
    assert((int) linearCutsNnz_ >= 0);
    assert(nele_cuts_jac == (int) linearCutsNnz_);
    if(iRow != NULL) {
      assert(jCol != NULL);
      assert(values == NULL);

      int offset = (index_style_ == Ipopt::TNLP::FORTRAN_STYLE);
	for(int i = 0; i < nele_cuts_jac; i++ ) {
	  iRow[i] = cutsiRow_[i] + m_ + offset; 
	  jCol[i] = cutsjCol_[i] + offset;
	}
    }
    else {
      assert(jCol == NULL);
      assert(values != NULL);
      IpBlasDcopy(nele_cuts_jac, cutsElems_, 1,
		  values,1);
    }
  }

  bool TMINLP2TNLP::eval_jac_g(Index n, const Number* x, bool new_x,
      Index m, Index nele_jac, Index* iRow,
      Index *jCol, Number* values)
  {
    bool return_code =
      tminlp_->eval_jac_g(n, x, new_x, m_, nele_jac - linearCutsNnz_,
			  iRow, jCol, values);
    if (return_code) {
      int start = nele_jac - linearCutsNnz_ ;
      if(iRow != NULL) iRow += start;
      if(jCol != NULL) jCol += start;
      if(values != NULL) values += start;
      eval_jac_g_add_linear_cuts(linearCutsNnz_ , iRow, jCol, 
                                 values);
    }
    return return_code;
  }

  bool TMINLP2TNLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess,
      Index* iRow, Index* jCol, Number* values)
  {
    return tminlp_->eval_h(n, x, new_x, obj_factor, m_ , lambda,
        new_lambda, nele_hess,
        iRow, jCol, values);
  }

  void TMINLP2TNLP::finalize_solution(SolverReturn status,
      Index n, const Number* x, const Number* z_L, const Number* z_U,
      Index m, const Number* g, const Number* lambda,
      Number obj_value,
      const IpoptData* ip_data,
      IpoptCalculatedQuantities* ip_cq)
  {
    assert(n==n_);
    assert(m==m_ + (int) nLinearCuts_);
    if (!x_sol_) {
      x_sol_ = new Number[n];
    }
    IpBlasDcopy(n, x, 1, x_sol_, 1);

    delete [] g_sol_;
    g_sol_ = new Number [m];
    IpBlasDcopy(m, g, 1, g_sol_, 1);
    delete [] duals_sol_;
    duals_sol_ = new Number[m + 2*n];
    if(lambda){
      IpBlasDcopy(m, lambda, 1, duals_sol_ + 2*n, 1);
      
      IpBlasDcopy(n, z_L, 1 , duals_sol_ , 1);
      IpBlasDcopy(n, z_U, 1 , duals_sol_ + n, 1);
    }

#if 0
    for (Index i=0; i<n; i++) {
      printf("x[%3d] = %13.6e z_L[%3d] = %13.6e z_U[%3d] = %13.6e\n", i, x[i], i, z_L[i], i, z_U[i]);
    }
    for (Index i=0; i<m; i++) {
      printf("lam[%3d] = %13.6e g[%3d] = %13.6e\n", i, lambda[i], i, g[i]);
    }
#endif

    return_status_ = status;
    obj_value_ = obj_value;

    if (IsValid(curr_warm_starter_)) {
      curr_warm_starter_->Finalize();
    }
  }

  void TMINLP2TNLP::throw_exception_on_bad_variable_bound(Index i)
  {
    assert(i >= 0 && i < n_);

    if (var_types_[i] == TMINLP::BINARY) {
      ASSERT_EXCEPTION(x_l_[i] == 0 && x_u_[i] == 1,
          TMINLP_INVALID_VARIABLE_BOUNDS,
          "Invalid variable bounds in TMINLP. All binaries must have 0,1 for variable bounds."
                      );
    }
    else if (var_types_[i] == TMINLP::INTEGER) {
      if(fabs(x_l_[i])<COIN_INT_MAX && fabs(x_u_[i]) < COIN_INT_MAX)//round to closest valid integer
      {
        int x_l = (int)ceil(x_l_[i]);
        int x_u = (int)floor(x_u_[i]);
        std::cerr<<"Inconsistent bounds on an integer"<<std::endl;
        ASSERT_EXCEPTION(x_l_[i] == (Number)x_l && x_u_[i] == (Number)x_u,
            TMINLP_INVALID_VARIABLE_BOUNDS,
            "Invalid variable bounds in TMINLP. All integer variables must have integer bounds."
                        );
      }
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
    const int &numcols = n_;
    const int &numrows = m_;

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
    for(int i=0 ; i < n_ ; i++) {
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
    //@{
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
    tminlp_->eval_upper_bound_f(n_, x, help);
    return help;
  }
  /** Resizes the starting point array in case cuts have been added or removed.
      Puts zeroes for eventually added cuts*/
  void 
  TMINLP2TNLP::resizeStartingPoint(int extra){
    //Always resize to have the capacity to hold twice the current number of constraints
    int newCapacity = 3 * n_ + 2 * (m_ + nLinearCuts_ + extra);
    double * new_x_init= new double [newCapacity];
    int end = min(capacity_x_init_, newCapacity);
    CoinCopyN(x_init_, end , new_x_init);
    int remaining = newCapacity - end;
    if(remaining){
      CoinZeroN(new_x_init + end , remaining);}
    delete [] x_init_;
    x_init_ = new_x_init;
    duals_init_ = x_init_ + n_;
    capacity_x_init_ = newCapacity;
  }

/** Add some cuts to the problem formulaiton (handles Quadratics).*/
void 
TMINLP2TNLP::addCuts(const OsiCuts &cuts){
   const OsiRowCut * first = cuts.rowCutPtr(0);
   addCuts(cuts.sizeRowCuts(), &first );
}
void
TMINLP2TNLP::addCuts(unsigned int numberCuts, const OsiRowCut ** cuts){

 int n,m,nnz_lag,nnz_hess;
 Ipopt::TNLP::IndexStyleEnum fort;
 get_nlp_info(n,m,nnz_lag,nnz_hess,fort);

 //count the number of non-zeroes
 //Cuts have to be added by rows for removeCuts to work
 unsigned int nnz = linearCutsNnz_;
 for(unsigned int i = 0 ; i < numberCuts ; i++)
 {
   nnz += cuts[i]->row().getNumElements();
 }
 int resizeNnz = nnz > linearCutsNnzCapacity_;
 int resizeCuts = nLinearCuts_ + numberCuts > linearCutsCapacity_;
 resizeLinearCuts(max((1+ resizeCuts) * linearCutsCapacity_, resizeCuts * numberCuts + linearCutsCapacity_), 
                  max((1 + resizeNnz) * linearCutsNnzCapacity_, resizeNnz * nnz + linearCutsNnzCapacity_));

 /* reinit nnz */
 nnz = linearCutsNnz_; 

 double * newCutsLowers = cutsLower_ + nLinearCuts_;
 double * newCutsUppers = cutsUpper_ + nLinearCuts_;

 for(unsigned int i = 0 ; i < numberCuts ; i++)
 {
   const int * ind = cuts[i]->row().getIndices();
   const double * values = cuts[i]->row().getElements();
   int size = cuts[i]->row().getNumElements();
   unsigned int iplusn = i + nLinearCuts_;
   for(int j = 0; j < size ; j++)
   {
    assert(ind[j] < n);
    cutsjCol_[nnz] = ind[j];
    cutsiRow_[nnz] = iplusn;
    cutsElems_[nnz++] = values[j];
   }
   newCutsLowers[i] = cuts[i]->lb();
   newCutsUppers[i] = cuts[i]->ub();
  }
 nLinearCuts_+=numberCuts;
 linearCutsNnz_=nnz;
}

void
TMINLP2TNLP::resizeLinearCuts(unsigned int newNumberCuts, unsigned int newNnz)
{
  if(newNumberCuts > linearCutsCapacity_)
  {
     double * newLower = new double[m_ + newNumberCuts];
     double * newUpper = new double[m_ + newNumberCuts];
     if(m_ + linearCutsCapacity_)
     {
       IpBlasDcopy(m_ + nLinearCuts_, g_l_, 1, newLower, 1);
       IpBlasDcopy(m_ + nLinearCuts_, g_u_, 1, newUpper, 1);
       delete [] g_l_;
       delete [] g_u_;
     }
     g_l_ = newLower;
     g_u_ = newUpper;
     cutsLower_ = newLower + m_;
     cutsUpper_ = newUpper + m_;
     linearCutsCapacity_ = newNumberCuts;
  }
  if(newNnz > linearCutsNnzCapacity_)
  {
    double * newElems = new double [newNnz];
    int * newiRow = new int [newNnz];
    int * newjCol = new int [newNnz];
    if(linearCutsNnzCapacity_)
    {
      IpBlasDcopy(linearCutsNnz_, cutsElems_, 1, newElems, 1);
      CoinCopyN(cutsiRow_, linearCutsNnz_, newiRow);
      CoinCopyN(cutsjCol_, linearCutsNnz_, newjCol);
      delete [] cutsElems_;
      delete [] cutsiRow_;
      delete [] cutsjCol_;
    }
    cutsElems_ = newElems;
    cutsiRow_ = newiRow;
    cutsjCol_ = newjCol;
    linearCutsNnzCapacity_ = newNnz;
  }
}
void
TMINLP2TNLP::removeCuts(unsigned int number, const int * toRemove){
   if(number==0) return;
    int n,m,nnz_lag,nnz_hess;
   Ipopt::TNLP::IndexStyleEnum fort;
   get_nlp_info(n,m,nnz_lag,nnz_hess,fort);

   int * sorted = new int[number];
   for(unsigned int i = 0 ; i < number ; i++) sorted[i] = toRemove[i] - m;
   std::sort(sorted, sorted + number);  
   int iNew = 0;
   int k = 0;
   for(unsigned int i = 0 ; i < linearCutsNnz_ ; i++)
   {
     if(sorted[k] < cutsiRow_[i]) k++;
     if(sorted[k] < cutsiRow_[i]) throw -1;
     if(sorted[k] == cutsiRow_[i]) continue;
     cutsiRow_[iNew] = cutsiRow_[i] - k;
     cutsjCol_[iNew] = cutsjCol_[i];
     cutsElems_[iNew++] = cutsElems_[i];
   }
 k=0;
 assert((int) nLinearCuts_ >= 0);
 for(int i = sorted[k] ; i < (int) nLinearCuts_ ; i++){
    if(sorted[k] < i) k++;
    if(sorted[k] < i) throw -1;
    if(sorted[k] == i) continue;
    cutsLower_[i - k] = cutsLower_[i];
    cutsUpper_[i - k] = cutsUpper_[i];
 }
 linearCutsNnz_ = iNew;
 nLinearCuts_ -= number;
 delete [] sorted;
}
void
TMINLP2TNLP::removeLastCuts(unsigned int number){
 number = nLinearCuts_ - number;
 assert((int) number >= 0);
 int iNew = 0;
   for(unsigned int i = 0 ; i < linearCutsNnz_ ; i++)
   {
     if(cutsiRow_[i] >= (int) number) {
     continue;}
     cutsiRow_[iNew] = cutsiRow_[i];
     cutsjCol_[iNew] = cutsjCol_[i];
     cutsElems_[iNew++] = cutsElems_[i];
   }
 linearCutsNnz_ = iNew;
 nLinearCuts_ = number;
 }


}// namespace Bonmin

