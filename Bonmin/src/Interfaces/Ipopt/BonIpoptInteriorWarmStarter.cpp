// (C) Copyright International Business Machines Corporation 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter          IBM    2006-03-09

#include "BonIpoptInteriorWarmStarter.hpp"
#include "IpDenseVector.hpp"

#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
namespace Bonmin
{
  IpoptInteriorWarmStarter::
  IpoptInteriorWarmStarter(Index n,
      const Number* x_l, const Number* x_u,
      Number nlp_lower_bound_inf,
      Number nlp_upper_bound_inf,
      bool store_several_iterates)
      :
      nlp_lower_bound_inf_(nlp_lower_bound_inf),
      nlp_upper_bound_inf_(nlp_upper_bound_inf),
      store_several_iterates_(store_several_iterates),
      n_(n),
      n_stored_iterates_(0)
  {
    x_l_prev_ = new double[n];
    x_u_prev_ = new double[n];

    for (Index i=0; i<n; i++) {
      x_l_prev_[i] = x_l[i];
      x_u_prev_[i] = x_u[i];
    }
  }

  IpoptInteriorWarmStarter::
  ~IpoptInteriorWarmStarter()
  {
    delete [] x_l_prev_;
    delete [] x_u_prev_;
  }

  bool IpoptInteriorWarmStarter::
  UpdateStoredIterates(AlgorithmMode mode,
      const IpoptData& ip_data,
      IpoptCalculatedQuantities& ip_cq)
  {
    // Don't store anything during the restoration phase
    if (mode==RestorationPhaseMode) {
      return true;
    }

    // Get some useful information out of the Ipopt objects
    Index iter = ip_data.iter_count();
    Number mu = ip_data.curr_mu();
    Number nlp_error = ip_cq.curr_nlp_error();
    Number primal_inf = ip_cq.curr_primal_infeasibility(NORM_MAX);
    Number dual_inf = ip_cq.curr_dual_infeasibility(NORM_MAX);
    Number complementarity = ip_cq.curr_complementarity(0., NORM_MAX);
    if (store_several_iterates_ || n_stored_iterates_==0) {
      // For now, we just store everything
      n_stored_iterates_++;
      stored_iter_.push_back(iter);
      stored_iterates_.push_back(ip_data.curr());
      stored_mu_.push_back(mu);
      stored_nlp_error_.push_back(nlp_error);
      stored_primal_inf_.push_back(primal_inf);
      stored_dual_inf_.push_back(dual_inf);
      stored_compl_.push_back(complementarity);
    }
    else {
      stored_iter_[0] = iter;
      stored_iterates_[0] = ip_data.curr();
      stored_mu_[0] = mu;
      stored_nlp_error_[0] = nlp_error;
      stored_primal_inf_[0] = primal_inf;
      stored_dual_inf_[0] = dual_inf;
      stored_compl_[0] = complementarity;
    }
    return true;
  }

  bool IpoptInteriorWarmStarter::
  Finalize()
  {
    // For now, nothing.  Later we could clean up, reduce storage etc.
    return true;
  }

  bool IpoptInteriorWarmStarter::
  WarmStartIterate(Index n, const Number* x_l_new,
      const Number* x_u_new,
      IteratesVector& warm_start_iterate)
  {
    assert(n==n_);

    if (n_stored_iterates_==0) {
      return false;
    }

    // For now let's just assume that we want to restore the 4-th latest
    // iterate from the previous solve

    Index iter_wanted = Max(0, n_stored_iterates_-5);

    SmartPtr<const Vector> prev_x = stored_iterates_[iter_wanted]->x();
    SmartPtr<const Vector> prev_s = stored_iterates_[iter_wanted]->s();
    SmartPtr<const Vector> prev_z_L = stored_iterates_[iter_wanted]->z_L();
    SmartPtr<const Vector> prev_z_U = stored_iterates_[iter_wanted]->z_U();
    SmartPtr<const Vector> prev_y_c = stored_iterates_[iter_wanted]->y_c();
    SmartPtr<const Vector> prev_y_d = stored_iterates_[iter_wanted]->y_d();
    SmartPtr<const Vector> prev_v_L = stored_iterates_[iter_wanted]->v_L();
    SmartPtr<const Vector> prev_v_U = stored_iterates_[iter_wanted]->v_U();

    const DenseVector* d_x = dynamic_cast<const DenseVector*> (GetRawPtr(prev_x));
    const DenseVector* d_s = dynamic_cast<const DenseVector*> (GetRawPtr(prev_s));
    const DenseVector* d_z_L = dynamic_cast<const DenseVector*> (GetRawPtr(prev_z_L));
    const DenseVector* d_z_U = dynamic_cast<const DenseVector*> (GetRawPtr(prev_z_U));
    const DenseVector* d_y_c = dynamic_cast<const DenseVector*> (GetRawPtr(prev_y_c));
    const DenseVector* d_y_d = dynamic_cast<const DenseVector*> (GetRawPtr(prev_y_d));
    const DenseVector* d_v_L = dynamic_cast<const DenseVector*> (GetRawPtr(prev_v_L));
    const DenseVector* d_v_U = dynamic_cast<const DenseVector*> (GetRawPtr(prev_v_U));

    const Number* x_vals_prev = d_x->Values();
    const Number* s_vals_prev = d_s->Values();
    const Number* z_L_vals_prev = d_z_L->Values();
    const Number* z_U_vals_prev = d_z_U->Values();
    const Number* y_c_vals_prev = d_y_c->Values();
    const Number* y_d_vals_prev = d_y_d->Values();
    const Number* v_L_vals_prev = d_v_L->Values();
    const Number* v_U_vals_prev = d_v_U->Values();

    DenseVector* d_x_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.x_NonConst()));
    DenseVector* d_s_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.s_NonConst()));
    DenseVector* d_z_L_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.z_L_NonConst()));
    DenseVector* d_z_U_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.z_U_NonConst()));
    DenseVector* d_y_c_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.y_c_NonConst()));
    DenseVector* d_y_d_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.y_d_NonConst()));
    DenseVector* d_v_L_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.v_L_NonConst()));
    DenseVector* d_v_U_new = dynamic_cast<DenseVector*> (GetRawPtr(warm_start_iterate.v_U_NonConst()));

    Number* x_vals_new = d_x_new->Values();
    Number* s_vals_new = d_s_new->Values();
    Number* z_L_vals_new = d_z_L_new->Values();
    Number* z_U_vals_new = d_z_U_new->Values();
    Number* y_c_vals_new = d_y_c_new->Values();
    Number* y_d_vals_new = d_y_d_new->Values();
    Number* v_L_vals_new = d_v_L_new->Values();
    Number* v_U_vals_new = d_v_U_new->Values();

    // Now copy the primal variables from the old to the new problem;
    // make sure that we take care of the fixed variables
    Index ix_prev = 0;
    Index ix_new = 0;
    Index izL_prev = 0;
    Index izL_new = 0;
    Index izU_prev = 0;
    Index izU_new = 0;
    for (Index i=0; i<n_; i++) {
      if (x_l_new[i]<x_u_new[i]) {
        DBG_ASSERT(x_l_prev_[i]<x_u_prev_[i]);
        x_vals_new[ix_new] = x_vals_prev[ix_prev];
        ix_new++;
        ix_prev++;
        if (x_l_new[i]>nlp_lower_bound_inf_) {
          DBG_ASSERT(x_l_prev_[i]>nlp_lower_bound_inf_);
          z_L_vals_new[izL_new] = z_L_vals_prev[izL_prev];
          izL_new++;
          izL_prev++;
        }
        if (x_u_new[i]<nlp_upper_bound_inf_) {
          DBG_ASSERT(x_u_prev_[i]<nlp_upper_bound_inf_);
          z_U_vals_new[izU_new] = z_U_vals_prev[izU_prev];
          izU_new++;
          izU_prev++;
        }
      }
      else if (x_l_prev_[i]<x_u_prev_[i]) {
        ix_prev++;
        izL_prev++;
        izU_prev++;
      }
    }
    DBG_ASSERT(ix_prev==prev_x->Dim());
    DBG_ASSERT(izL_prev==prev_z_L->Dim());
    DBG_ASSERT(izU_prev==prev_z_U->Dim());
    DBG_ASSERT(ix_new==warm_start_iterate.x()->Dim());
    DBG_ASSERT(izL_new==warm_start_iterate.z_L()->Dim());
    DBG_ASSERT(izU_new==warm_start_iterate.z_U()->Dim());

    // Now copy all the values for the iterates that don't change in dimension
    DBG_ASSERT(prev_s->Dim()==warm_start_iterate.s()->Dim());
    DBG_ASSERT(prev_y_d->Dim()==warm_start_iterate.s()->Dim());
    DBG_ASSERT(prev_y_d->Dim()==warm_start_iterate.y_d()->Dim());
    for (Index i=0; i<prev_s->Dim(); i++) {
      s_vals_new[i] = s_vals_prev[i];
      y_d_vals_new[i] = y_d_vals_prev[i];
    }
    DBG_ASSERT(prev_y_c->Dim()==warm_start_iterate.y_c()->Dim());
    for (Index i=0; i<prev_y_c->Dim(); i++) {
      y_c_vals_new[i] = y_c_vals_prev[i];
    }
    DBG_ASSERT(prev_v_L->Dim()==warm_start_iterate.v_L()->Dim());
    for (Index i=0; i<prev_v_L->Dim(); i++) {
      v_L_vals_new[i] = v_L_vals_prev[i];
    }
    DBG_ASSERT(prev_v_U->Dim()==warm_start_iterate.v_U()->Dim());
    for (Index i=0; i<prev_v_U->Dim(); i++) {
      v_U_vals_new[i] = v_U_vals_prev[i];
    }

    return true;
  }
}

