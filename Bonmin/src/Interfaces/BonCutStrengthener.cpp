// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:   Andreas Waechter                 IBM    2007-03-30

#include "BonCutStrengthener.hpp"
#include "IpBlas.hpp"

namespace Bonmin
{
  using namespace Ipopt;

  CutStrengthener::CutStrengthener(SmartPtr<TNLPSolver> tnlp_solver)
    :
    tnlp_solver_(tnlp_solver)
  {
    tnlp_solver_->Options()->clear();
    tnlp_solver_->Initialize("strength.opt");
    tnlp_solver_->Options()->SetStringValue("hessian_approximation","limited-memory");
    tnlp_solver_->Options()->SetStringValue("mu_strategy", "adaptive");
  }

  CutStrengthener::~CutStrengthener()
  {
  }

  bool CutStrengthener::StrengthenCut(SmartPtr<TMINLP> tminlp,
				      int constr_index,
				      CoinPackedVector& row,
				      int n,
				      const double* x,
				      const double* x_l,
				      const double* x_u,
				      double& lb,
				      double& ub)
  {
    // Get data to set up the NLP to be solved
    Index nele_grad_gi;
    Index* jCol = new Index[n];
    bool new_x = true;
    if (!tminlp->eval_grad_gi(n, x, new_x, constr_index, nele_grad_gi,
				jCol, NULL)) {
      delete [] jCol;
      return false;
    }

    bool lower_bound;
    if (lb <= -DBL_MAX) {
      assert(ub < DBL_MAX);
      lower_bound = false;
    }
    else {
      assert(ub >= DBL_MAX);
      lower_bound = true;
    }
    SmartPtr<StrengtheningTNLP> stnlp =
      new StrengtheningTNLP(tminlp, row, lower_bound, n, x, x_l, x_u,
			    constr_index, nele_grad_gi, jCol);

    delete [] jCol;

    TNLPSolver::ReturnStatus status =
      tnlp_solver_->OptimizeTNLP(GetRawPtr(stnlp));

    if (status == TNLPSolver::solvedOptimal ||
	status == TNLPSolver::solvedOptimalTol) {
      const Number tiny_move = 0e-8;
      const Number final_bound = stnlp->StrengthenedBound();
      if (lower_bound) {
	lb = final_bound - tiny_move;
      }
      else {
	ub = final_bound + tiny_move;
      }
    }
    else {
      return false;
    }

    return true;
  }

  CutStrengthener::StrengtheningTNLP::
  StrengtheningTNLP(SmartPtr<TMINLP> tminlp,
		    const CoinPackedVector& cut,
		    bool lower_bound,
		    Index n,
		    const Number* starting_point,
		    const double* x_l_orig,
		    const double* x_u_orig,
		    Index constr_index,
		    Index nvar_constr /** Number of variables in constraint */,
		    const Index* jCol)
    :
    tminlp_(tminlp),
    n_orig_(n),
    constr_index_(constr_index),
    nvar_constr_(nvar_constr),
    lower_bound_(lower_bound),
    have_final_bound_(false)
  {
    starting_point_ = new Number[n_orig_];
    x_full_ = new Number[n_orig_];
    IpBlasDcopy(n_orig_, starting_point, 1, starting_point_, 1);
    IpBlasDcopy(n_orig_, starting_point, 1, x_full_, 1);

    obj_grad_ = new Number[nvar_constr_];
    x_l_ = new Number[nvar_constr_];
    x_u_ = new Number[nvar_constr_];
    const Number zero = 0.;
    IpBlasDcopy(nvar_constr_, &zero, 0, obj_grad_, 1);

    const int cut_nele = cut.getNumElements();
    const int* cut_indices = cut.getIndices();
    const double* cut_elements = cut.getElements();

    for (int i = 0; i<cut_nele; i++) {
      const int& idx = cut_indices[i];

      // ToDo: This could be done more efficiently
      Index jidx = -1;
      for (int j=0; j<nvar_constr_; j++) {
	if (idx == jCol[j]) {
	  jidx = j;
	  break;
	}
      }
      if (jidx < 0) {
	std::cerr << "There is an index in the cut that does not appear in the constraint.\n";
	exit(-99);
      }

      if (lower_bound) {
	obj_grad_[jidx] = cut_elements[i];
      }
      else {
	obj_grad_[jidx] = -cut_elements[i];
      }

    }

    var_indices_ = new Index[nvar_constr_];
    for (int i=0; i<nvar_constr_; i++) {
      const Index& j = jCol[i];
      var_indices_[i] = j;
      x_l_[i] = x_l_orig[j];
      x_u_[i] = x_u_orig[j];
    }
  }

  CutStrengthener::StrengtheningTNLP::~StrengtheningTNLP()
  {
    delete [] obj_grad_;
    delete [] x_l_;
    delete [] x_u_;
    delete [] var_indices_;
    delete [] starting_point_;
    delete [] x_full_;
  }

  bool CutStrengthener::StrengtheningTNLP::
  get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	       Index& nnz_h_lag, IndexStyleEnum& index_style)
  {
    n = nvar_constr_;
    m = 1;
    nnz_jac_g = nvar_constr_;
    nnz_h_lag = 0;
    index_style = C_STYLE;

    Index n_orig;
    Index nnz_jac_g_orig;
    Index nnz_h_lag_orig;
    TNLP::IndexStyleEnum index_style_orig;
    if(!tminlp_->get_nlp_info(n_orig, m_orig_, nnz_jac_g_orig, nnz_h_lag_orig,
			      index_style_orig)) {
      return false;
    }
    if (n_orig_ != n_orig) {
      std::cerr << "Number of variables inconsistent in StrengtheningTNLP::get_nlp_info\n";
      return false;
    }

    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  get_bounds_info(Index n, Number* x_l, Number* x_u,
		  Index m, Number* g_l, Number* g_u)
  {
    Number* x_l_orig = new Number[n_orig_];
    Number* x_u_orig = new Number[n_orig_];
    Number* g_l_orig = new Number[m_orig_];
    Number* g_u_orig = new Number[m_orig_];

    if (!tminlp_->get_bounds_info(n_orig_, x_l_orig, x_u_orig,
				  m_orig_, g_l_orig, g_u_orig)) {
      return false;
    }

    g_l[0] = g_l_orig[constr_index_];
    g_u[0] = g_u_orig[constr_index_];

    for (Index i=0; i<nvar_constr_; i++) {
      x_l[i] = x_l_[i];
      x_u[i] = x_u_[i];
    }

    delete [] x_l_orig;
    delete [] x_u_orig;
    delete [] g_l_orig;
    delete [] g_u_orig;

    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  get_starting_point(Index n, bool init_x, Number* x,
		     bool init_z, Number* z_L, Number* z_U,
		     Index m, bool init_lambda,
		     Number* lambda)
  {
    assert(!init_z && !init_lambda);
    assert(n = nvar_constr_);
    if (init_x) {
      for (Index i=0; i<n; i++) {
	x[i] = starting_point_[var_indices_[i]];
      }
    }

    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    obj_value = 0.;
    for (Index i=0; i<n; i++) {
      obj_value += obj_grad_[i]*x[i];
    }
    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    IpBlasDcopy(n, obj_grad_, 1, grad_f, 1);
    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    update_x_full(x);
    return tminlp_->eval_gi(n_orig_, x_full_, new_x, constr_index_, g[0]);
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_jac_g(Index n, const Number* x, bool new_x,
	     Index m, Index nele_jac, Index* iRow, Index *jCol,
	     Number* values)
  {
    if (iRow) {
      DBG_ASSERT(jCol && !values);
      for (Index i=0; i<nele_jac; i++) {
	iRow[i] = 0;
	jCol[i] = i;
      }
    }
    else {
      DBG_ASSERT(!iRow && values);
      update_x_full(x);
      if (!tminlp_->eval_grad_gi(n_orig_, x_full_, new_x, constr_index_,
				 nele_jac, NULL, values)) {
	return false;
      }
    }
    return true;
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_h(Index n, const Number* x, bool new_x,
	 Number obj_factor, Index m, const Number* lambda,
	 bool new_lambda, Index nele_hess,
	 Index* iRow, Index* jCol, Number* values)
  {
    std::cerr << "At the moment, StrengtheningTNLP::eval_h is not yet implemented\n";
    return false;
  }

  void CutStrengthener::StrengtheningTNLP::
  update_x_full(const Number *x)
  {
    for (Index i=0; i<nvar_constr_; i++) {
      x_full_[var_indices_[i]] = x[i];
    }
  }

  void CutStrengthener::StrengtheningTNLP::
  finalize_solution(SolverReturn status,
		    Index n, const Number* x, const Number* z_L, const Number* z_U,
		    Index m, const Number* g, const Number* lambda,
		    Number obj_value,
		    const IpoptData* ip_data,
		    IpoptCalculatedQuantities* ip_cq)
  {
    if (status == SUCCESS || status == STOP_AT_ACCEPTABLE_POINT) {
      have_final_bound_ = true;
      strengthened_bound_ = obj_value;
    }
    else {
      have_final_bound_ = false;
    }
  }

  Number CutStrengthener::StrengtheningTNLP::StrengthenedBound() const
  {
    DBG_ASSERT(have_final_bound_);
    if (lower_bound_) {
      return strengthened_bound_;
    }
    else {
      return -strengthened_bound_;
    }
  }
  

} // namespace Bonmin
