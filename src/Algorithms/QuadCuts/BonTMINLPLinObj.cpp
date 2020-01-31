// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 08/16/2007

#include "BonTMINLPLinObj.hpp"

using namespace Ipopt;

namespace Bonmin{
   /** Default constructor*/
   TMINLPLinObj::TMINLPLinObj():
   tminlp_(NULL), m_(0), n_(0), 
   nnz_jac_(0)
   {}

   void 
   TMINLPLinObj::gutsOfDestructor(){
     tminlp_ = NULL;
   }
 
   
   TMINLPLinObj::~TMINLPLinObj(){
   gutsOfDestructor();}


  void
   TMINLPLinObj::setTminlp(SmartPtr<TMINLP> tminlp){
      gutsOfDestructor();
      tminlp_ = tminlp;
      int n,m, nnz_jac, nnz_h;
      Ipopt::TNLP::IndexStyleEnum index_style;
      tminlp_->get_nlp_info(n, m , nnz_jac, nnz_h, index_style);
      offset_ = index_style == Ipopt::TNLP::FORTRAN_STYLE;
      n_ = n+1;
      m_ = m+1;
      nnz_jac_ = nnz_jac + n + 1;
  }

   bool 
   TMINLPLinObj::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                         Index& nnz_h_lag, 
                                         Ipopt::TNLP::IndexStyleEnum& index_style){
     assert(IsValid(tminlp_));

     bool return_value = 
          tminlp_->get_nlp_info(n,m,nnz_jac_g, nnz_h_lag,index_style);
     m = m_;
     n = n_;
     nnz_jac_g = nnz_jac_;
     return return_value;
   }
   

   bool 
   TMINLPLinObj::get_scaling_parameters(Number& obj_scaling,
                                                   bool& use_x_scaling, Index n,
                                                   Number* x_scaling,
                                                   bool& use_g_scaling, Index m,
                                                   Number* g_scaling){
     assert(IsValid(tminlp_));
     assert(m == m_);


     if(g_scaling && use_g_scaling)
       g_scaling[0] = 1.;
     if(x_scaling && use_x_scaling)
       x_scaling[n - 1] = 1.;
     obj_scaling = 1.;
     double dummy = 1.;
     double * mod_obj_scaling = &dummy;
     if(use_g_scaling && g_scaling){
       mod_obj_scaling = g_scaling;}
     return tminlp_->get_scaling_parameters(*mod_obj_scaling, use_x_scaling, n - 1, x_scaling,
                                                         use_g_scaling, m - 1, g_scaling + 1);
   }





   bool
   TMINLPLinObj::get_constraints_linearity(Index m, 
					   Ipopt::TNLP::LinearityType* const_types){
      assert(IsValid(tminlp_));
      assert(m == m_ );
              const_types[0] = Ipopt::TNLP::NON_LINEAR;
        return tminlp_->get_constraints_linearity(m-1, const_types +1);}

   bool
   TMINLPLinObj::get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u){
     assert(IsValid(tminlp_));
     assert(m == m_);
     assert(n == n_);

     x_l[n-1] = -DBL_MAX;
     x_u[n-1] = DBL_MAX;
    
     g_l[0] = -DBL_MAX;
     g_u[0] = 0.;
     return tminlp_->get_bounds_info(n - 1, x_l, x_u, m_ - 1,
                                     g_l +1, g_u + 1);
   }

   bool
   TMINLPLinObj::get_starting_point(Index n, bool init_x, Number* x,
                                                bool init_z, Number* z_L, Number* z_U,
                                                Index m, bool init_lambda,
                                                Number* lambda){
     assert(IsValid(tminlp_));
     assert(m == m_);

     bool return_value = tminlp_->get_starting_point(n - 1, init_x, x, init_z, z_L, z_U,
                                            m - 1, init_lambda, lambda + 1);
     tminlp_->eval_f(n-1, x, true, x[n-1]);
     if(init_lambda && lambda != NULL)
       lambda[0] = 0;
     return return_value;
   }
   

   bool
   TMINLPLinObj::eval_g(Index n, const Number* x, bool new_x,
                                   Index m, Number* g){
     assert(IsValid(tminlp_));
     assert(m == m_);
     assert(n == n_);

     bool ret_val =  tminlp_->eval_f(n - 1, x, new_x, g[0]);
     g[0] -= x[n -1];
     return ret_val && tminlp_->eval_g(n - 1, x, false, m - 1, g+1);
   }

   bool
   TMINLPLinObj::eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values){
     assert(IsValid(tminlp_));
     assert(m == m_);
     assert(n == n_);
     assert(nele_jac == nnz_jac_);
     bool ret_val = true;
     if(values == NULL)
     {
       for(int i = 0 ; i < n_ ; i++){
         iRow[i] = offset_; jCol[i] = i + offset_;}
       bool ret_val = tminlp_->eval_jac_g(n -1, x, new_x, m_ -1, nnz_jac_ - n_, iRow + n_, jCol + n_, NULL);
       for(int i = n_ ; i < nnz_jac_ ; i++){//shift by 1
         iRow[i]++;}
       return ret_val;
     }
     else {
       ret_val &= tminlp_->eval_grad_f(n-1, x, new_x, values);
       values[n-1] = -1;
       ret_val &= tminlp_->eval_jac_g(n - 1, x, false, m - 1, nele_jac - n_, NULL, NULL, values + n);
     }
     return ret_val;
   }

   bool
   TMINLPLinObj::eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values){
     assert(IsValid(tminlp_));
     assert(m == m_);
     assert(n == n_);
     return tminlp_->eval_h(n_ - 1, x, new_x, (lambda != NULL)? lambda[0]: 1., m_ - 1, (lambda != NULL)? lambda + 1: NULL, new_lambda,
                                  nele_hess, iRow, jCol, values);
    }


   bool
   TMINLPLinObj:: eval_gi(Index n, const Number* x, bool new_x,
			 Index i, Number& gi){
     assert(IsValid(tminlp_));
     assert(i < m_);
     assert(n == n_);
     if(i == 0){
      bool ret_val = tminlp_->eval_f(n-1, x, new_x, gi);
      gi -= x[n -1];
      return ret_val;}
      else
      return tminlp_->eval_gi(n-1, x, new_x, i - 1, gi);
   }

   bool
   TMINLPLinObj::eval_grad_gi(Index n, const Number* x, bool new_x,
			      Index i, Index& nele_grad_gi, Index* jCol,
			      Number* values){
     assert(IsValid(tminlp_));
     assert(i < m_);
     assert(n == n_);
     if(i == 0){
      if(jCol != NULL){
        for(int i = 0 ; i < n ; i++) jCol[i] = i + offset_;
      }
      bool ret_val= tminlp_->eval_grad_f(n-1, x, new_x, values);
      values[n-1] = -1;
      return ret_val;}
      else
        return tminlp_->eval_grad_gi(n - 1, x, new_x, i - 1, nele_grad_gi,
                        jCol, values);
   }


}/* Ends namespace Bonmin.*/
