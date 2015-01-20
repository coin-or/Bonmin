// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 05/25/2005


#include "BonTNLP2FPNLP.hpp"

using namespace Ipopt;

namespace Bonmin
{
  TNLP2FPNLP::TNLP2FPNLP(const SmartPtr<TNLP> tnlp, double objectiveScalingFactor):
      tnlp_(tnlp),
      inds_(),
      vals_(),
      lambda_(1.),
      sigma_(1.),
      norm_(2),
      objectiveScalingFactor_(objectiveScalingFactor),
      use_feasibility_pump_objective_(false),
      use_cutoff_constraint_(false),
      use_local_branching_constraint_(false),
      cutoff_(COIN_DBL_MAX),
      rhs_local_branching_constraint_(COIN_DBL_MAX),
      index_style_(TNLP::C_STYLE)
  {}

  TNLP2FPNLP::TNLP2FPNLP(const SmartPtr<TNLP> tnlp, const SmartPtr<TNLP2FPNLP> other):
      tnlp_(tnlp),
      inds_(other->inds_),
      vals_(other->vals_),
      lambda_(other->lambda_),
      sigma_(other->sigma_),
      norm_(other->norm_),
      objectiveScalingFactor_(other->objectiveScalingFactor_),
      use_feasibility_pump_objective_(other->use_feasibility_pump_objective_),
      use_cutoff_constraint_(other->use_cutoff_constraint_),
      use_local_branching_constraint_(other->use_local_branching_constraint_),
      cutoff_(other->cutoff_),
      rhs_local_branching_constraint_(other->rhs_local_branching_constraint_),
      index_style_(other->index_style_)
  {}

  TNLP2FPNLP::~TNLP2FPNLP()
  {
  }

  void
  TNLP2FPNLP::set_cutoff(Number cutoff)
  {
    Number epsilon = 1.0e-6;
    if(cutoff > 1.0e-8)
      cutoff_ = (1-epsilon) * cutoff;
    else if(cutoff < -1.0e-8)
      cutoff_ = (1+epsilon) * cutoff;
    else
      cutoff_ = - epsilon;
  }

  void
  TNLP2FPNLP::set_dist_to_point_obj(size_t n, const Number * vals, const Index * inds)
  {
    inds_.resize(n);
    vals_.resize(n);
    std::copy(vals, vals + n, vals_.begin());
    std::copy(inds, inds + n, inds_.begin());
  }

  /** Compute the distance to the current point to which distance is minimized. */
  double
  TNLP2FPNLP::dist_to_point(const Number *x)
  {
    double ret_val = 0;
    assert(vals_.size() == inds_.size());
    if(norm_ == 2){
      for(unsigned int i = 0; i < vals_.size() ; i++) {
        ret_val += ( x[inds_[i]] - vals_[i] ) * ( x[inds_[i]] - vals_[i] );
      }
    }
    else if(norm_ == 1){
      for(unsigned int i = 0 ; i < vals_.size() ; i++) {
        if(vals_[i] <= 0.1)
         ret_val += x[inds_[i]];
        else{
	  ret_val += (1.0 - x[inds_[i]]);
	 //         ret_val -= x[inds_[i]];
        }
      }
    }
    return ret_val;
  }

  bool
  TNLP2FPNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
      Index& nnz_h_lag,
      TNLP::IndexStyleEnum& index_style)
  {
    bool ret_code = tnlp_->get_nlp_info(n, m , nnz_jac_g, nnz_h_lag,
        index_style);

    index_style_ = index_style; // this function is called before any
                                // other that uses index_style_

    if(use_feasibility_pump_objective_ && norm_ == 2)
      nnz_h_lag += (int)vals_.size();

    if(use_cutoff_constraint_ && use_local_branching_constraint_) {
      m += 2;
      nnz_jac_g += (n + (int)vals_.size());
    }
    else if(use_cutoff_constraint_) {
      m++;
      nnz_jac_g += n;
    }
    else if(use_local_branching_constraint_) {
      m++;
      nnz_jac_g += (int)vals_.size();
    }

    return ret_code;
  }

  bool 
  TNLP2FPNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
			      Index m, Number* g_l, Number* g_u)
  {
    bool ret_code;

    if(use_cutoff_constraint_ && use_local_branching_constraint_) {
      ret_code = tnlp_->get_bounds_info(n, x_l , x_u, m-2, g_l, g_u);
      g_l[m-2] = - COIN_DBL_MAX;
      g_u[m-2] = cutoff_;
      g_l[m-1] = - COIN_DBL_MAX;
      g_u[m-1] = rhs_local_branching_constraint_;
    }
    else if(use_cutoff_constraint_) {
      ret_code = tnlp_->get_bounds_info(n, x_l , x_u, m-1, g_l, g_u);
      g_l[m-1] = - COIN_DBL_MAX;
      g_u[m-1] = cutoff_;
    }
    else if(use_local_branching_constraint_) {
      ret_code = tnlp_->get_bounds_info(n, x_l , x_u, m-1, g_l, g_u);
      g_l[m-1] = - COIN_DBL_MAX;
      g_u[m-1] = rhs_local_branching_constraint_;
    }
    else
      ret_code = tnlp_->get_bounds_info(n, x_l , x_u, m, g_l, g_u);

    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_f(Index n, const Number* x, bool new_x,
      Number& obj_value)
  {
    bool ret_code = tnlp_->eval_f(n, x, new_x, obj_value);//for new_x

    if(use_feasibility_pump_objective_) {
      obj_value *= (1 - lambda_) * sigma_;
      obj_value += objectiveScalingFactor_*lambda_*dist_to_point(x);
    }

    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_grad_f(Index n, const Number* x, bool new_x,
      Number* grad_f)
  {
    bool ret_code = tnlp_->eval_grad_f(n, x, new_x, grad_f);

    if(use_feasibility_pump_objective_) {
      for(int i = 0 ; i < n ; i++) {
	grad_f[i] *= (1-lambda_) * sigma_;
      }
      if(norm_ ==2){
	for(unsigned int i = 0 ; i < inds_.size() ; i++) {
	  grad_f[inds_[i]] += objectiveScalingFactor_*2*lambda_*( x[inds_[i]] - vals_[i] );
	}
      }
      else{
	for(unsigned int i = 0 ; i < inds_.size() ; i++) {
	  if(vals_[i] <= 0.1)
	    grad_f[inds_[i]] += objectiveScalingFactor_*lambda_;
	  else
	    grad_f[inds_[i]] -= objectiveScalingFactor_*lambda_;
	}
      }
    }
   
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_g(Index n, const Number* x, bool new_x,
		     Index m, Number* g)
  {
    bool ret_code;

    if(use_cutoff_constraint_ && use_local_branching_constraint_) {
      ret_code = tnlp_->eval_g(n, x, new_x, m-2, g);
      // compute the value of the cutoff constraint
      Number obj_value;
      if(eval_f(n, x, new_x, obj_value))
	g[m-2] = obj_value;
      else
	ret_code = false;
      // compute the value of the local branching constraint
      Number g_local_branching = 0.0;
      for(unsigned int i = 0 ; i < vals_.size() ; i++) {
        if(vals_[i] <= 0.1)
          g_local_branching += x[inds_[i]];
        else
	  g_local_branching += (1.0 - x[inds_[i]]);
      }
      g[m-1] = g_local_branching;
    }
    else if(use_cutoff_constraint_) {
      ret_code = tnlp_->eval_g(n, x, new_x, m-1, g);
      Number obj_value;
      if(eval_f(n, x, new_x, obj_value))
	g[m-1] = obj_value;
      else
	ret_code = false;
    }
    else if(use_local_branching_constraint_) {
      ret_code = tnlp_->eval_g(n, x, new_x, m-1, g);
      Number g_local_branching = 0.0;
      for(unsigned int i = 0 ; i < vals_.size() ; i++) {
        if(vals_[i] <= 0.1)
          g_local_branching += x[inds_[i]];
        else
	  g_local_branching += (1.0 - x[inds_[i]]);
      }
      g[m-1] = g_local_branching;
    }
    else
      ret_code = tnlp_->eval_g(n, x, new_x, m, g);
    
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_jac_g(Index n, const Number* x, bool new_x,
			 Index m, Index nele_jac, Index* iRow,
			 Index *jCol, Number* values)
  {
    bool ret_code;

    if(use_cutoff_constraint_ && use_local_branching_constraint_) {
      int n_integers = (int)vals_.size();
      ret_code = tnlp_->eval_jac_g(n, x, new_x, m, nele_jac - n - n_integers, 
				   iRow, jCol, values);

      if (iRow && jCol && !values) { //Initialization phase 
	int index_correction = (index_style_ == TNLP::C_STYLE) ? 0 : 1;
	// compute the jacobian contribution of the cutoff constraint
	int k = nele_jac - n - n_integers;
	iRow += k;
	jCol += k;
	for(int i = 0; i< n; i++) {
	  iRow[i] = m - 2 + index_correction;
	  jCol[i] = i + index_correction;
	}
	// compute the jacobian contribution of the local branching constraint
	k = nele_jac - n_integers;
	iRow += k;
	jCol += k;
	for(int i = 0; i< n_integers; i++) {
	  iRow[i] = m - 1 + index_correction;
	  jCol[i] = inds_[i] + index_correction;
	}
      }
      else if (!iRow & !jCol && values) { //computation phase
	// compute the jacobian contribution of the cutoff constraint
	Number* grad_f = new Number[n];
	bool ret_code_grad_f = eval_grad_f(n, x, new_x, grad_f);
	if(ret_code_grad_f) {
	  int k = nele_jac - n - n_integers;
	  values += k;
	  for(int i = 0; i< n; i++) {
	    values[i] = grad_f[i];
	  }
	}
	else
	  ret_code = false;
	delete [] grad_f;
	// compute the jacobian contribution of the local branching constraint
	int k = nele_jac - n_integers;
	values += k;
	for(int i = 0; i< n_integers; i++) {
	  if(vals_[i] <= 0.1)
	    values[i] = 1.0;
	  else
	    values[i] = -1.0;
	}
      }	  
      else { //error phase
	DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
      }
    }
    else if(use_cutoff_constraint_) {
      ret_code = tnlp_->eval_jac_g(n, x, new_x, m, nele_jac - n, 
				   iRow, jCol, values);
      
      if (iRow && jCol && !values) { //Initialization phase 
	int index_correction = (index_style_ == TNLP::C_STYLE) ? 0 : 1;
	int k = nele_jac - n;
	iRow += k;
	jCol += k;
	for(int i = 0; i< n; i++) {
	  iRow[i] = m - 1 + index_correction;
	  jCol[i] = i + index_correction;
	}
      }
      else if (!iRow & !jCol && values) { //computation phase
	Number* grad_f = new Number[n];
	bool ret_code_grad_f = eval_grad_f(n, x, new_x, grad_f);
	if(ret_code_grad_f) {
	  int k = nele_jac - n;
	  values += k;
	  for(int i = 0; i< n; i++) {
	    values[i] = grad_f[i];
	  }
	}
	else
	  ret_code = false;
	delete [] grad_f;
      }	  
      else { //error phase
	DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
      }
    }
    else if(use_local_branching_constraint_) {
      int n_integers = (int)vals_.size();
      ret_code = tnlp_->eval_jac_g(n, x, new_x, m, nele_jac - n_integers, 
				   iRow, jCol, values);
      
      if (iRow && jCol && !values) { //Initialization phase 
	int index_correction = (index_style_ == TNLP::C_STYLE) ? 0 : 1;
	int k = nele_jac - n_integers;
	iRow += k;
	jCol += k;
	for(int i = 0; i< n_integers; i++) {
	  iRow[i] = m - 1 + index_correction;
	  jCol[i] = inds_[i] + index_correction;
	}
      }
      else if (!iRow & !jCol && values) { //computation phase
	int k = nele_jac - n_integers;
	values += k;
	for(int i = 0; i< n_integers; i++) {
	  if(vals_[i] <= 0.1)
	    values[i] = 1.0;
	  else
	    values[i] = -1.0;
	}
      }	  
      else { //error phase
	DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
      }
    }
    else
      ret_code = tnlp_->eval_jac_g(n, x, new_x, m, nele_jac, 
				   iRow, jCol, values);
    
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_h(Index n, const Number* x, bool new_x,
		     Number obj_factor, Index m, const Number* lambda,
		     bool new_lambda, Index nele_hess,
		     Index* iRow, Index* jCol, Number* values)
  {
    bool ret_code;

    int  nnz_obj_h = (norm_ == 2) ? (int)inds_.size() : 0;

    if(use_cutoff_constraint_ && use_local_branching_constraint_) {
      double coef_obj = (iRow != NULL)?0 : lambda[m - 2];
      ret_code = tnlp_->eval_h(n, x, new_x, obj_factor*(1-lambda_)*sigma_ + coef_obj, 
			       m - 2, lambda, new_lambda, nele_hess - nnz_obj_h, 
			       iRow, jCol, values);
    }
    else if(use_cutoff_constraint_) {
      double coef_obj = (iRow != NULL)?0 : lambda[m - 1];
      ret_code = tnlp_->eval_h(n, x, new_x, obj_factor*(1-lambda_)*sigma_ + coef_obj, 
			       m - 1, lambda, new_lambda, nele_hess - nnz_obj_h, 
			       iRow, jCol, values);
    }
    else if(use_local_branching_constraint_) {
      ret_code = tnlp_->eval_h(n, x, new_x, obj_factor*(1-lambda_)*sigma_, 
			       m - 1, lambda, new_lambda, nele_hess - nnz_obj_h, 
			       iRow, jCol, values);
    }
    else { // this is the original feasibility pump implementation
      ret_code = tnlp_->eval_h(n, x, new_x, obj_factor*(1-lambda_)*sigma_, 
			       m, lambda, new_lambda, nele_hess - nnz_obj_h, 
			       iRow, jCol, values);
    }

    //Now add extra elements corresponding to the hessian of the distance
    if(use_feasibility_pump_objective_ && norm_ == 2){
      if (iRow && jCol && !values) //Initialization phase
	{
    int index_correction = (index_style_ == TNLP::C_STYLE) ? 0 : 1;
	  int k = nele_hess - nnz_obj_h;
	  iRow += k;
	  jCol += k;
	  for(unsigned int i = 0; i < inds_.size() ; i++)
	    {
	      iRow[i] = inds_[i] + index_correction;
	      jCol[i] = inds_[i] + index_correction;
	    }
	  DBG_ASSERT(k==nele_hess);
	}
      else if (!iRow & !jCol && values) //computation phase
	{
	  int k = nele_hess - nnz_obj_h;
	  values += k;
	  for(unsigned int i = 0; i < inds_.size() ; i++)
	    {
	      values[i] = 2* objectiveScalingFactor_* lambda_* obj_factor;
	    }
	  DBG_ASSERT(k==nele_hess);
	}
      else //error phase
	{
	  DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
	}
    }

    return ret_code;
  }

  void
  TNLP2FPNLP::finalize_solution(SolverReturn status,
      Index n, const Number* x, const Number* z_L, const Number* z_U,
      Index m, const Number* g, const Number* lambda,
      Number obj_value,
      const IpoptData* ip_data,
      IpoptCalculatedQuantities* ip_cq)
  {
      int m2 = m;
      if(use_cutoff_constraint_) {
        m2--;
      }
      if(use_local_branching_constraint_) {
        m2--;
      }
    tnlp_->finalize_solution(status,n, x, z_L, z_U,m2, g, lambda, obj_value,
			     ip_data, ip_cq);
  }
}

