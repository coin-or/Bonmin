// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 05/25/2005


#include "BonTNLP2FPNLP.hpp"
#include "IpBlas.hpp"
namespace Bonmin
{
  TNLP2FPNLP::TNLP2FPNLP(const SmartPtr<TNLP> tnlp, double objectiveScalingFactor):
      tnlp_(tnlp),
      inds_(),
      vals_(),
      lambda_(1.),
      sigma_(1.),
      norm_(2),
      objectiveScalingFactor_(objectiveScalingFactor)
  {}

  TNLP2FPNLP::~TNLP2FPNLP()
  {
  }


  void
  TNLP2FPNLP::set_dist2point_obj(int n, const Number * vals, const Index * inds)
  {
    inds_.resize(n);
    vals_.resize(n);
    IpBlasDcopy(n, vals, 1, vals_(), 1);
    CoinCopyN(inds, n, inds_());
  }

  /** Compute the distance to the current point to which distance is minimized. */
  double
  TNLP2FPNLP::dist2point(const Number *x)
  {
    double ret_val = 0;
    assert(norm_ > 0 && norm < 3);
    assert(vals_.size() == inds_.size());
    if(norm_ == 2){
      for(unsigned int i = 0; i < vals_.size() ; i++) {
        ret_val += ( x[inds_[i]] - vals_[i] ) * ( x[inds_[i]] - vals_[i] );
      }
    }
    else if(norm_ == 1){
      for(unsigned int i = 0 ; i < vals_.size() ; i++) {
        ret_val += fabs(x[inds_[i]] - vals_[i]);
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
    if(norm_ == 2)
      nnz_h_lag += vals_.size();
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_f(Index n, const Number* x, bool new_x,
      Number& obj_value)
  {
    bool ret_code = tnlp_->eval_f(n, x, new_x, obj_value);//for new_x
    obj_value *= (1 - lambda_) * sigma_;
    obj_value = objectiveScalingFactor_*lambda_*dist2point(x);
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_grad_f(Index n, const Number* x, bool new_x,
      Number* grad_f)
  {
    bool ret_code = tnlp_->eval_grad_f(n, x, new_x, grad_f);
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
        if(x[inds_[i]] - vals_[i] >= 0.)
          grad_f[inds_[i]] += objectiveScalingFactor_*lambda_;
        else
          grad_f[inds_[i]] -= objectiveScalingFactor_*lambda_;
      }
    }
    
   
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess,
      Index* iRow, Index* jCol, Number* values)
  {

    int  nnz_obj_h = inds_.size();
    //Call the function for the Hessian of the original constraint system
    bool ret_code = tnlp_->eval_h(n, x, new_x, obj_factor*(1-lambda_)*sigma_, 
                                  m, lambda, new_lambda, nele_hess - nnz_obj_h, 
                                  iRow, jCol, values);
    //Now add extra elements corresponding to the hessian of the distance
    if(norm_ == 2){
      if (iRow && jCol && !values) //Initialization phase
      {
        int k = nele_hess - nnz_obj_h;
        iRow += k;
        jCol += k;
        for(unsigned int i = 0; i < inds_.size() ; i++)
        {
          iRow[i] = inds_[i] + 1;
          jCol[i] = inds_[i] + 1;
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
    tnlp_->finalize_solution(status,n, x, z_L, z_U,m, g, lambda, obj_value,
			     ip_data, ip_cq);
  }
}

