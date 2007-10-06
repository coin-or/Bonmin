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
      var_types_(NULL),
      x_l_(NULL),
      x_u_(NULL),
      orig_x_l_(NULL),
      orig_x_u_(NULL),
      g_l_(NULL),
      g_u_(NULL),
      x_init_(NULL),
      capacity_x_init_(0),
      x_init_user_(NULL),
      x_sol_(NULL),
      g_sol_(NULL),
      duals_sol_(NULL),
//      return_status_(NOT_SOLVED),
      obj_value_(1e100),
      curr_warm_starter_(NULL),
      need_new_warm_starter_(true)
  {
    // read the nlp size and bounds information from
    // the TMINLP and keep an internal copy. This way the
    // caller can modify the bounds that are sent to Ipopt;
    DBG_ASSERT(IsValid(tminlp_));

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
    duals_init_(NULL),
    capacity_x_init_(other.capacity_x_init_),
    x_sol_(NULL),
    g_sol_(NULL),
    duals_sol_(NULL)
  {
    var_types_ = new TMINLP::VariableType[n_];
    for (Index i=0; i<n_; i++) {
      var_types_[i] = other.var_types_[i];
    }

    x_l_ = new Number[n_];
    x_u_ = new Number[n_]; // Those are copied in copyUserModification

    orig_x_l_ = new Number[n_];
    orig_x_u_ = new Number[n_];
    IpBlasDcopy(n_, other.orig_x_l_, 1, orig_x_l_, 1);
    IpBlasDcopy(n_, other.orig_x_u_, 1, orig_x_u_, 1);

    g_l_ = new Number[m_];
    g_u_ = new Number[m_];
    IpBlasDcopy(m_, other.g_l_, 1, g_l_, 1);
    IpBlasDcopy(m_, other.g_u_, 1, g_u_, 1);

    x_init_ = new Number[capacity_x_init_];
    x_init_user_ = new Number[n_];
    IpBlasDcopy(n_, other.x_init_user_, 1, x_init_user_, 1);

    copyUserModification(other);
  }

  TMINLP2TNLP::~TMINLP2TNLP()
  {
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
  }

  /** Copies the modification made to TNLP by the user (modifications
      such as changing bound changing starting point,...).
      this and other should be two instances of the same problem
      I am trying to mimic a copy construction for Cbc
      use with great care not safe.
  */
  void
  TMINLP2TNLP::copyUserModification(const TMINLP2TNLP& other)
  {
    DBG_ASSERT(x_l_);
    DBG_ASSERT(x_u_);
    DBG_ASSERT(other.x_l_);
    DBG_ASSERT(other.x_u_);
    DBG_ASSERT(n_ == other.n_);
    DBG_ASSERT(m_ == other.m_);

    IpBlasDcopy(n_, other.x_l_, 1, x_l_, 1);
    IpBlasDcopy(n_, other.x_u_, 1, x_u_, 1);

    if(other.duals_init_) {
      duals_init_ = &x_init_[n_];
      IpBlasDcopy(capacity_x_init_, other.x_init_, 1, x_init_, 1);
    }
    else
      IpBlasDcopy(n_, other.x_init_, 1, x_init_, 1);
    return_status_ = other.return_status_;

    if(other.x_sol_ !=NULL) {
      //DBG_ASSERT(return_status_ != NOT_SOLVED);
      Set_x_sol(n_,other.x_sol_);
    }

    if(other.g_sol_!=NULL) {
//	DBG_ASSERT(return_status_ != NOT_SOLVED);
      g_sol_ = new Number [m_];
      IpBlasDcopy(m_, other.g_sol_, 1, g_sol_, 1);
    }

    if(other.duals_sol_!=NULL) {
//	DBG_ASSERT(return_status_ != NOT_SOLVED);
      duals_sol_ = new Number[m_ + 2*n_];
      IpBlasDcopy(2*n_+ m_, other.duals_sol_, 1, duals_sol_, 1);
    }

    obj_value_ = other.obj_value_;

    curr_warm_starter_ = other.curr_warm_starter_;

    nlp_lower_bound_inf_ = other.nlp_lower_bound_inf_;

    nlp_upper_bound_inf_ = other.nlp_upper_bound_inf_;

    need_new_warm_starter_ = other.need_new_warm_starter_;
  }

  void TMINLP2TNLP::SetVariablesBounds(Index n,
                                       const Number * x_l,
                                       const Number * x_u)
  {
    DBG_ASSERT(n==n_);
    IpBlasDcopy(n_, x_l, 1, x_l_, 1);
    IpBlasDcopy(n_, x_u, 1, x_u_, 1);
  }

   void TMINLP2TNLP::SetVariablesLowerBounds(Index n,
                                       const Number * x_l)
  {
    DBG_ASSERT(n==n_);
    IpBlasDcopy(n_, x_l, 1, x_l_, 1);
  }

   void TMINLP2TNLP::SetVariablesUpperBounds(Index n,
                                       const Number * x_u)
  {
    DBG_ASSERT(n==n_);
    IpBlasDcopy(n_, x_u, 1, x_u_, 1);
  }

  void TMINLP2TNLP::SetVariableBounds(Index var_no, Number x_l, Number x_u)
  {
    DBG_ASSERT(var_no >= 0 && var_no < n_);
    x_l_[var_no] = x_l;
    x_u_[var_no] = x_u;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetVariableLowerBound(Index var_no, Number x_l)
  {
    DBG_ASSERT(var_no >= 0 && var_no < n_);
    x_l_[var_no] = x_l;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetVariableUpperBound(Index var_no, Number x_u)
  {
    DBG_ASSERT(var_no >= 0 && var_no < n_);
    x_u_[var_no] = x_u;
    //    throw_exception_on_bad_variable_bound(var_no);
  }

  void TMINLP2TNLP::SetStartingPoint(Index n,const Number* x_init)
  {
    DBG_ASSERT(n == n_);
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
    DBG_ASSERT(n == n_);
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
    DBG_ASSERT(m == m_ + 2*n_ );

    if(!duals_init_)
      duals_init_ = x_init_ + n_;

    IpBlasDcopy(m, duals_init, 1, duals_init_, 1);

  }

  /** Set the contiuous solution */
  void TMINLP2TNLP::Set_x_sol(Index n, const Number* x_sol)
  {
    DBG_ASSERT(n == n_);
    if (!x_sol_) {
      x_sol_ = new Number[n];
    }
    IpBlasDcopy(n_, x_sol, 1, x_sol_, 1);
  }

  /** Change the type of the variable */
  void TMINLP2TNLP::SetVariableType(Index n, TMINLP::VariableType type)
  {
    DBG_ASSERT(n >= 0 && n < n_);
    var_types_[n] = type;
  }

  bool TMINLP2TNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
      Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    n = n_;
    m = m_ + tminlp_->nLinearCuts_;
    nnz_jac_g = nnz_jac_g_ + tminlp_->linearCutsNnz_;
    nnz_h_lag = nnz_h_lag_;
    index_style = index_style_;
    return true;
  }

  bool TMINLP2TNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
      Index m, Number* g_l, Number* g_u)
  {
    DBG_ASSERT(n==n_);
    DBG_ASSERT(m==m_);
    IpBlasDcopy(n_, x_l_, 1, x_l, 1);
    IpBlasDcopy(n_, x_u_, 1, x_u, 1);
    IpBlasDcopy(m_, g_l_, 1, g_l, 1);
    IpBlasDcopy(m_, g_u_, 1, g_u, 1);
    IpBlasDcopy(tminlp_->nLinearCuts_, tminlp_->cutsLower_, 1, &g_l[m_],1); 
    IpBlasDcopy(tminlp_->nLinearCuts_, tminlp_->cutsUpper_, 1, &g_u[m_],1); 
    return true;
  }

  bool TMINLP2TNLP::get_starting_point(Index n, bool init_x, Number* x,
      bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda,
      Number* lambda)
  {
    assert(m==m_ + tminlp_->nLinearCuts_);
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
      int size= 3*n + m_ + tminlp_->nLinearCuts_;
      if(size > capacity_x_init_)
        resizeStartingPoint();
      IpBlasDcopy(m_ + tminlp_->nLinearCuts_ , duals_init_ + 2*n , 1, lambda, 1);
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
    DBG_ASSERT(nCuts == nCuts);
    DBG_ASSERT(n == n_);
    // Add the linear cuts
    CoinZeroN(g , tminlp_->nLinearCuts_);
    for(int nnz = 0 ; nnz < tminlp_->linearCutsNnz_ ; nnz++){
      g[tminlp_->cutsiRow_[nnz]] += tminlp_->cutsElems_[nnz] * x[tminlp_->cutsjCol_[nnz]];
    }

  }

  bool TMINLP2TNLP::eval_g(Index n, const Number* x, bool new_x,
      Index m, Number* g)
  {
    int return_code = tminlp_->eval_g(n, x, new_x, m_, g);
    if (return_code) {
      eval_g_add_linear_cuts(n, x, tminlp_->nLinearCuts_, 
                             g + m_);
    }
    return return_code;
  }

  void TMINLP2TNLP::eval_jac_g_add_linear_cuts(Index nele_cuts_jac, Index* iRow,
					       Index* jCol, Number* values)
  {
    DBG_ASSERT(nele_cuts_jac == tminlp_->linearCutsNnz_);
    if(iRow != NULL) {
      DBG_ASSERT(jCol != NULL);
      DBG_ASSERT(values == NULL);

      int offset = (index_style_ == TNLP::FORTRAN_STYLE);
	for(int i = 0; i < nele_cuts_jac; i++ ) {
	  iRow[i] = tminlp_->cutsiRow_[i] + m_ + offset; 
	  jCol[i] = tminlp_->cutsjCol_[i] + offset;
	}
    }
    else {
      DBG_ASSERT(jCol == NULL);
      DBG_ASSERT(values != NULL);
      IpBlasDcopy(nele_cuts_jac, tminlp_->cutsElems_, 1,
		  values,1);
    }
  }

  bool TMINLP2TNLP::eval_jac_g(Index n, const Number* x, bool new_x,
      Index m, Index nele_jac, Index* iRow,
      Index *jCol, Number* values)
  {
    bool return_code =
      tminlp_->eval_jac_g(n, x, new_x, m_, nele_jac - tminlp_->linearCutsNnz_,
			  iRow, jCol, values);
    if (return_code) {
      int start = nele_jac - tminlp_->linearCutsNnz_ ;
      if(iRow != NULL) iRow += start;
      if(jCol != NULL) jCol += start;
      if(values != NULL) values += start;
      eval_jac_g_add_linear_cuts(tminlp_->linearCutsNnz_ , iRow, jCol, 
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
    assert(m==m_ + tminlp_->nLinearCuts_);
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
    DBG_ASSERT(i >= 0 && i < n_);

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
  void 
  TMINLP2TNLP::evaluateUpperBoundingFunction(const double * x){
    tminlp_->eval_upper_bound_f(n_, x, obj_value_);
  }
  /** Resizes the starting point array in case cuts have been added or removed.
      Puts zeroes for eventually added cuts*/
  void 
  TMINLP2TNLP::resizeStartingPoint(){
    //Always resize to have the capacity to hold twice the current number of constraints
    int newCapacity = 3 * n_ + 2 * (m_ + tminlp_->nLinearCuts_);
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
}// namespace Ipopt

