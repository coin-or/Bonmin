// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2008
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                   (derived from BonTMINLP2TNLP.cpp)            12/22/2006
// Authors :


#include "BonBranchingTQP.hpp"
#include "IpBlas.hpp"
#include "IpAlgTypes.hpp"
#include <string>
#include <fstream>
#include <sstream>

using namespace Ipopt;

namespace Bonmin
{
  BranchingTQP::BranchingTQP(SmartPtr<TMINLP2TNLP> tminlp2tnlp)
    :
    tminlp2tnlp_(tminlp2tnlp)
  {
    bool retval = tminlp2tnlp_->get_nlp_info(n_, m_, nnz_jac_g_,
					     nnz_h_lag_, index_style_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't get NLP infor in BranchingTQP");
    //DBG_ASSERT(x_sol_);
    //DBG_ASSERT(duals_sol_);

    obj_grad_ = new Number[n_];
    obj_hess_ = new Number[nnz_h_lag_];
    obj_hess_irow_ = new Index[nnz_h_lag_];
    obj_hess_jcol_ = new Index[nnz_h_lag_];
    g_vals_ = new Number[m_];
    g_jac_ = new Number[nnz_jac_g_];
    g_jac_irow_ = new Index[nnz_jac_g_];
    g_jac_jcol_ = new Index[nnz_jac_g_];

    const Number* x_sol = tminlp2tnlp_->x_sol();
    const Number* duals_sol = tminlp2tnlp_->duals_sol();

    // Compute all nonlinear values at the starting point so that we
    // have all the information for the QP
    bool new_x = true;   // ToDo: maybe NOT new?
    retval = tminlp2tnlp_->eval_f(n_, x_sol, new_x, obj_val_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate objective function in BranchingTQP");
    new_x = false;
    retval = tminlp2tnlp_->eval_grad_f(n_, x_sol, new_x, obj_grad_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate objective gradient in BranchingTQP");
    bool new_lambda = true; // ToDo: maybe NOT new?
    retval = tminlp2tnlp_->eval_h(n_, x_sol, new_x, 1., m_, duals_sol + 2 * n_,
			     new_lambda, nnz_h_lag_, obj_hess_irow_,
			     obj_hess_jcol_, NULL);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate objective Hessian structure in BranchingTQP");
    if (index_style_ == TNLP::FORTRAN_STYLE) {
      for (Index i=0; i<nnz_h_lag_; i++) {
	obj_hess_irow_[i]--;
	obj_hess_jcol_[i]--;
      }
    }
    retval = tminlp2tnlp_->eval_h(n_, x_sol, new_x, 1., m_, duals_sol + 2*n_,
			     new_lambda, nnz_h_lag_, NULL, NULL, obj_hess_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate objective Hessian values in BranchingTQP");
    retval = tminlp2tnlp_->eval_g(n_, x_sol, new_x, m_, g_vals_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate constraint values in BranchingTQP");
    retval = tminlp2tnlp_->eval_jac_g(n_, x_sol, new_x, m_, nnz_jac_g_,
				 g_jac_irow_, g_jac_jcol_, NULL);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate constraint Jacobian structure in BranchingTQP");
    if (index_style_ == TNLP::FORTRAN_STYLE) {
      for (Index i=0; i<nnz_jac_g_; i++) {
	g_jac_irow_[i]--;
	g_jac_jcol_[i]--;
      }
    }
    retval = tminlp2tnlp_->eval_jac_g(n_, x_sol, new_x, m_, nnz_jac_g_,
				 NULL, NULL, g_jac_);
    ASSERT_EXCEPTION(retval, TMINLP_INVALID,
		     "Can't evaluate constraint Jacobian values in BranchingTQP");

    // Keep copy of original x_sol and duals_sol values
    x_sol_copy_ = new Number[n_];
    IpBlasDcopy(n_, x_sol, 1, x_sol_copy_, 1);
    duals_sol_copy_ = new Number[m_ + 2*n_];
    IpBlasDcopy(m_+2*n_, duals_sol, 1, duals_sol_copy_, 1);
#if 0
    for (int i=0; i<n_; i++) {
      printf("x_sol_copy_[%3d] = %15.8e duals_sol_copy_[%3d] = %15.8e obj_grad_[%3d] = %15.8e\n", i,x_sol_copy_[i],i,duals_sol_copy_[i],i,obj_grad_[i]);
    }
    for (int i=0; i<m_; i++) {
      printf("g_vals_[%3d] = %15.8e\n", i, g_vals_[i]);
    }
    for (int i=0; i<nnz_h_lag_; i++) {
      printf("Hess[%3d,%3d] = %15.8e\n",obj_hess_irow_[i],obj_hess_jcol_[i],obj_hess_[i]);
    }
    for (int i=0; i<nnz_jac_g_; i++) {
      printf("Jac[%3d,%3d] = %15.8e\n",g_jac_irow_[i],g_jac_jcol_[i],g_jac_[i]);
    }
#endif
  }

  BranchingTQP::~BranchingTQP()
  {
    delete [] obj_grad_;
    delete [] obj_hess_;
    delete [] obj_hess_irow_;
    delete [] obj_hess_jcol_;
    delete [] g_vals_;
    delete [] g_jac_;
    delete [] g_jac_irow_;
    delete [] g_jac_jcol_;
    delete [] x_sol_copy_;
    delete [] duals_sol_copy_;
  }

  bool BranchingTQP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				  Index& nnz_h_lag,
				  IndexStyleEnum& index_style)
  {
    n = n_;
    m = m_;
    nnz_jac_g = nnz_jac_g_;
    nnz_h_lag = nnz_h_lag_;
    index_style = index_style_;
    return true;
  }

  bool BranchingTQP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u)
  {
    DBG_ASSERT(n == n_);
    DBG_ASSERT(m == m_);
    bool retval = tminlp2tnlp_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
    // correct for displacement
    for (int i=0; i<n; i++) {
      x_l[i] -= x_sol_copy_[i];
      x_u[i] -= x_sol_copy_[i];
    }
    // include the right hand side of the constraint
    for (int i=0; i<m; i++) {
      g_l[i] -= g_vals_[i];
      g_u[i] -= g_vals_[i];
    }
    return retval;
  }

  bool BranchingTQP::get_starting_point(Index n, bool init_x, Number* x,
      bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda,
      Number* lambda)
  {
    DBG_ASSERT(n==n_);
    if (init_x == true) {
      const double zero = 0.;
      IpBlasDcopy(n, &zero, 0, x, 1);
    }
    if (init_z == true) {
      DBG_ASSERT(duals_sol_copy_);
      IpBlasDcopy(n, duals_sol_copy_, 1, z_L, 1);
      IpBlasDcopy(n, duals_sol_copy_ + n, 1, z_U, 1);

    }
    if(init_lambda == true) {
      DBG_ASSERT(duals_sol_copy_);
      IpBlasDcopy(m_, duals_sol_copy_ + 2*n_, 1, lambda, 1);
      for(int i = m_ ; i < m; i++)
      {
        lambda[i] = 0.;
      }
    }

    return true;
  }

  bool
  BranchingTQP::get_constraints_linearity(Index m, LinearityType* const_types)
  {
    for (int i=0; i<m; i++) {
      const_types[i] = LINEAR;
    }
    return true;
  }

  bool BranchingTQP::eval_f(Index n, const Number* x, bool new_x,
      Number& obj_value)
  {
    DBG_ASSERT(n == n_);

    obj_value = IpBlasDdot(n, x, 1, obj_grad_, 1);
    for (int i=0; i<nnz_h_lag_; i++) {
      Index& irow = obj_hess_irow_[i];
      Index& jcol = obj_hess_jcol_[i];
      if (irow!=jcol) {
	obj_value += obj_hess_[i]*x[irow]*x[jcol];
      }
      else {
	obj_value += 0.5*obj_hess_[i]*x[irow]*x[irow];
      }
    }

    return true;
  }

  bool BranchingTQP::eval_grad_f(Index n, const Number* x, bool new_x,
				 Number* grad_f)
  {
    DBG_ASSERT(n == n_);

    IpBlasDcopy(n_, obj_grad_, 1, grad_f, 1);
    for (int i=0; i<nnz_h_lag_; i++) {
      Index& irow = obj_hess_irow_[i];
      Index& jcol = obj_hess_jcol_[i];
      grad_f[irow] += obj_hess_[i]*x[jcol];
      if (irow!=jcol) {
	grad_f[jcol] += obj_hess_[i]*x[irow];
      }
    }

    return true;
  }

  bool BranchingTQP::eval_g(Index n, const Number* x, bool new_x,
			    Index m, Number* g)
  {
    DBG_ASSERT(n == n_);

    const double zero = 0.;
    IpBlasDcopy(m_, &zero, 0, g, 1);
    for (Index i=0; i<nnz_jac_g_; i++) {
      Index& irow = g_jac_irow_[i];
      Index& jcol = g_jac_jcol_[i];
      g[irow] += g_jac_[i]*x[jcol];
    }

    return true;
  }

  bool BranchingTQP::eval_jac_g(Index n, const Number* x, bool new_x,
				Index m, Index nele_jac, Index* iRow,
				Index *jCol, Number* values)
  {
    if (iRow != NULL) {
      DBG_ASSERT(jCol != NULL);
      DBG_ASSERT(values == NULL);
      if (index_style_ == TNLP::FORTRAN_STYLE) {
	for (Index i=0; i<nnz_jac_g_; i++) {
	  iRow[i] = g_jac_irow_[i] + 1;
	  jCol[i] = g_jac_jcol_[i] + 1;
	}
      }
      else {
	for (Index i=0; i<nnz_jac_g_; i++) {
	  iRow[i] = g_jac_irow_[i];
	  jCol[i] = g_jac_jcol_[i];
	}
      }
    }
    else {
      IpBlasDcopy(nnz_jac_g_, g_jac_, 1, values, 1);
    }

    return true;
  }

  bool BranchingTQP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess,
      Index* iRow, Index* jCol, Number* values)
  {
    DBG_ASSERT(nele_hess == nnz_h_lag_);

    if (iRow != NULL) {
      DBG_ASSERT(jCol != NULL);
      DBG_ASSERT(values == NULL);
      if (index_style_ == TNLP::FORTRAN_STYLE) {
	for (Index i=0; i<nele_hess; i++) {
	  iRow[i] = obj_hess_irow_[i] + 1;
	  jCol[i] = obj_hess_jcol_[i] + 1;
	}
      }
      else {
	for (Index i=0; i<nele_hess; i++) {
	  iRow[i] = obj_hess_irow_[i];
	  jCol[i] = obj_hess_jcol_[i];
	}
      }
    }
    else {
      if (obj_factor==0.) {
	const Number zero = 0.;
	IpBlasDcopy(nele_hess, &zero, 0, values, 1);
      }
      else {
	IpBlasDcopy(nele_hess, obj_hess_, 1, values, 1);
	if (obj_factor != 1.) {
	  IpBlasDscal(nele_hess, obj_factor, values, 1);
	}
      }
    }

    return true;
  }

  void BranchingTQP::finalize_solution(SolverReturn status,
				       Index n, const Number* x,
				       const Number* z_L, const Number* z_U,
				       Index m, const Number* g,
				       const Number* lambda, Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq)
  {
    // Translate displacement solution back to a solution in the real
    // variables
    double* xx = new double[n];
    for (int i=0; i<n; i++) {
      xx[i] = x_sol_copy_[i] + x[i];
    }

    double obj = obj_value + obj_val_;
    if(status == Ipopt::LOCAL_INFEASIBILITY) 
    obj = obj_value;
    tminlp2tnlp_->finalize_solution(status, n, xx, z_L, z_U, m, g, lambda,
				    obj, ip_data, ip_cq);
    delete [] xx;
  }

}
// namespace Ipopt

