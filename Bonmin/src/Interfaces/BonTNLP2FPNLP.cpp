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
      n_(-1),
      nMax_(-1),
      inds_(NULL),
      vals_(NULL),
      objectiveScalingFactor_(objectiveScalingFactor)
  {}

  TNLP2FPNLP::~TNLP2FPNLP()
  {
    if(inds_!=NULL)
      delete [] inds_;
    inds_=NULL;
    if(vals_!=NULL)
      delete [] vals_;
    vals_=NULL;
    nMax_ = -1;

  }


  void
  TNLP2FPNLP::set_dist2point_obj(int n, const Number * vals, const Index * inds)
  {
    if(n > nMax_)//resize tables
    {
      if(inds_!=NULL)
        delete [] inds_;
      if(vals_!=NULL)
        delete [] vals_;
      inds_ = new int[n];
      vals_ = new double[n];
      nMax_ = n;
    }
    //copy data
    n_ = n;
    IpBlasDcopy(n_, vals, 1, vals_, 1);
    for(int i = 0 ; i < n_ ; i++)
      inds_[i] = inds[i];
  }

  /** Compute the norm-2 distance to the current point to which distance is minimized. */
  double
  TNLP2FPNLP::dist2point(const Number *x)
  {
    double ret_val = 0;
    for(int i = 0; i < n_ ; i++) {
      ret_val += ( x[inds_[i]] - vals_[i] ) * ( x[inds_[i]] - vals_[i] );
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
    nnz_h_lag += n_;
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_f(Index n, const Number* x, bool new_x,
      Number& obj_value)
  {
    bool ret_code = tnlp_->eval_f(n, x, new_x, obj_value);//for new_x
    obj_value = objectiveScalingFactor_*dist2point(x);
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_grad_f(Index n, const Number* x, bool new_x,
      Number* grad_f)
  {
    int n1,m,d1,d2;
    TNLP::IndexStyleEnum index_style;
    get_nlp_info(n1,m,d1,d2,index_style);
    assert(n==n1);
    bool ret_code = tnlp_->eval_grad_f(n, x, new_x, grad_f);
    for(int i = 0 ; i < n ; i++) {
      grad_f[i] = 0.;
    }
    for(int i = 0 ; i < n_ ; i++) {
      grad_f[inds_[i]] = objectiveScalingFactor_*2 *( x[inds_[i]] - vals_[i] );
    }
    return ret_code;
  }

  bool
  TNLP2FPNLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess,
      Index* iRow, Index* jCol, Number* values)
  {

    int  nnz_obj_h = n_;
    //Call the function for the Hessian of the original constraint system
    bool ret_code = tnlp_->eval_h(n, x, new_x, 0., m, lambda, new_lambda, nele_hess - nnz_obj_h, iRow, jCol, values);
    //Now add extra elements corresponding to the hessian of the distance
    if (iRow && jCol && !values) //Initialization phase
    {
      int k = nele_hess - nnz_obj_h;
      for(int i = 0; i < n_ ; i++)
      {
        iRow[k] = inds_[i] + 1;
        jCol[k] = inds_[i] + 1;
        k++;
      }
      DBG_ASSERT(k==nele_hess);
    }
    else if (!iRow & !jCol && values) //computation phase
    {
      int k = nele_hess - nnz_obj_h;
      for(int i = 0; i < n_ ; i++)
      {
        values[k] = 2* objectiveScalingFactor_* obj_factor;
        k++;
      }
      DBG_ASSERT(k==nele_hess);
    }
    else //error phase
    {
      DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
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

