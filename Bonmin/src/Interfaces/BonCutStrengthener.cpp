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

  CutStrengthener::CutStrengthener(SmartPtr<TNLPSolver> tnlp_solver,
				   SmartPtr<OptionsList> options)
    :
    tnlp_solver_(tnlp_solver)
  {
    options->GetIntegerValue("oa_log_level", oa_log_level_, tnlp_solver->prefix());
    options->GetEnumValue("cut_strengthening_type", cut_strengthening_type_,
			  tnlp_solver->prefix());
    options->GetEnumValue("disjunctive_cut_type", disjunctive_cut_type_,
			  tnlp_solver->prefix());

    tnlp_solver_->options()->clear();
    if (!tnlp_solver_->Initialize("strength.opt")) {
      throw CoinError("CutStrengthener","CutStrengthener","Error during initialization of tnlp_solver_");
    }
    tnlp_solver_->options()->SetStringValue("hessian_approximation","limited-memory");
    tnlp_solver_->options()->SetStringValue("mu_strategy", "adaptive");
  }

  CutStrengthener::~CutStrengthener()
  {
  }

  bool CutStrengthener::HandleOneCut(bool is_tight, TMINLP* tminlp,
				     TMINLP2TNLP* problem,
				     const double* minlp_lb,
				     const double* minlp_ub,
				     const int gindex, CoinPackedVector& cut,
				     double& cut_lb, double& cut_ub,
				     int n, const double* x,
				     double infty)
  {
    bool retval = true;
    const int cut_nele = cut.getNumElements();
    const int* cut_indices = cut.getIndices();
    const TMINLP::VariableType* vartypes = problem->var_types();
    const double* cut_elements = cut.getElements();
    switch (disjunctive_cut_type_) {
    case DC_None:
      if (!is_tight) {
	retval = StrengthenCut(tminlp, gindex, cut, n, x, minlp_lb, minlp_ub,
			       cut_lb, cut_ub);
      }
      break;
    case DC_MostFractional: {
      // First check if there is a discrete variable to do disjunction on
      int imostfra = -1;
      double viol = 1e-6;
      for (int i=0; i<cut_nele; i++) {
	const int& idx = cut_indices[i];
	if (idx < n && (vartypes[idx] == TMINLP::BINARY ||
			vartypes[idx] == TMINLP::INTEGER)) {
	  const double& xi = x[idx];
	  const double this_viol = CoinMin(xi-floor(xi),ceil(xi)-xi);
	  if (this_viol > viol) {
	    imostfra = i;
	    viol = this_viol;
	  }
	}
      }
      if (imostfra == -1) {
	// No disjunction to be done
	if (!is_tight) {
	  retval = StrengthenCut(tminlp, gindex, cut, n, x, minlp_lb, minlp_ub,
				 cut_lb, cut_ub);
	}
      }
      else {
	//Do the disjunction for imostfra
	const int& idx = cut_indices[imostfra];
	const double& xi = x[idx];
	if (oa_log_level_>=2) {
	  printf("Doing disjunction for constr %d on x[%d] = %e\n", gindex, idx, xi);
	}
	const double down_xi = floor(xi);
	double* changed_bnds = new double[n];
	// First the down disjunction:
	CoinCopyN(minlp_ub, n, changed_bnds);
	changed_bnds[idx] = down_xi;
	double cut_lb_down = cut_lb;
	double cut_ub_down = cut_ub;
	retval = StrengthenCut(tminlp, gindex, cut, n, x, minlp_lb,
			       changed_bnds, cut_lb_down, cut_ub_down);
	double cut_lb_up = cut_lb;
	double cut_ub_up = cut_ub;
	if (retval) {
	  CoinCopyN(minlp_lb, n, changed_bnds);
	  changed_bnds[idx] = down_xi + 1.;
	  retval = StrengthenCut(tminlp, gindex, cut, n, x, changed_bnds,
				 minlp_ub, cut_lb_up, cut_ub_up);
	}
	delete [] changed_bnds;
	if (retval) {
	  // change the cut
	  double coeff_xi = cut_elements[imostfra];
	  const double old_coeff = coeff_xi;
	  if (cut_lb <= -infty) {
	    // cut has upper bound
	    coeff_xi += cut_ub_down - cut_ub_up;
	    cut_ub = cut_ub_up + (down_xi+1.)*(cut_ub_down - cut_ub_up);
	  }
	  else {
	    // cut has lower bound
	    coeff_xi += cut_lb_down - cut_lb_up;
	    cut_lb = cut_lb_up + (down_xi+1.)*(cut_lb_down - cut_lb_up);
	  }
	  cut.setElement(imostfra, coeff_xi);
	  printf("old coeff = %e new = %e\n", old_coeff, coeff_xi);
	}
      }
    }
      break;
    default:
      std::cerr << "Invalid case for disjunctive_cut_type_ in CutStrengthener HandleOneCut\n";
      exit(-2);
      break;
    }

    return retval;
  }

  bool CutStrengthener::ComputeCuts(OsiCuts &cs,
				    TMINLP* tminlp,
				    TMINLP2TNLP* problem,
				    const int gindex, CoinPackedVector& cut,
				    double& cut_lb, double& cut_ub,
				    const double g_val, const double g_lb,
				    const double g_ub,
				    int n, const double* x,
				    double infty)
  {
    //printf("before: lb = %e ub = %e rl = %e ru = %e g = %e\n", lb[i], ub[i], rowLower[bindi], rowUpper[bindi], g[bindi]);
    // First check if the cut is indeed away from the constraint
    bool is_tight = false;
    if (gindex==-1) {
      is_tight = true;
    }
    else {
      const Number tight_tol = 1e-8;
      if (cut_lb <= -infty && g_ub - g_val <= tight_tol) {
	is_tight = true;
      }
      else if (cut_ub >= infty && g_val - g_lb <= tight_tol) {
	is_tight = true;
      }
    }
    if (cut_strengthening_type_ == CS_StrengthenedGlobal ||
	cut_strengthening_type_ == CS_StrengthenedGlobal_StrengthenedLocal) {
      const double orig_lb = cut_lb;
      const double orig_ub = cut_ub;
      bool retval = HandleOneCut(is_tight, tminlp, problem,
				 problem->orig_x_l(),
				 problem->orig_x_u(), gindex, cut,
				 cut_lb, cut_ub, n, x, infty);
      if (!retval) {
        if (oa_log_level_ >= 1) {
          printf(" Error during strengthening of global cut for constraint %d\n", gindex);
        }
      }
      else {
        if (oa_log_level_ >=2 && (fabs(orig_lb-cut_lb)>1e-4 ||
				  fabs(orig_ub-cut_ub)>1e-4)) {
	  if (orig_ub < infty) {
	    printf(" Strengthening ub of global cut for constraint %d from %e to %e\n", gindex, orig_ub, cut_ub);
	  }
	  else {
	    printf(" Strengthening lb of global cut for constraint %d from %e to %e\n", gindex, orig_lb, cut_lb);
	  }
        }
      }
    }
    if (cut_strengthening_type_ == CS_UnstrengthenedGlobal_StrengthenedLocal ||
	cut_strengthening_type_ == CS_StrengthenedGlobal_StrengthenedLocal) {
      Number lb2 = cut_lb;
      Number ub2 = cut_ub;
      CoinPackedVector cut2(cut);
      bool retval = HandleOneCut(is_tight, tminlp, problem, problem->x_l(),
				 problem->x_u(), gindex, cut2,
				 lb2, ub2, n, x, infty);
      if (!retval) {
        if (oa_log_level_ >= 1) {
          printf(" Error during strengthening of local cut for constraint %d\n", gindex);
        }
      }
      else {
        const Number localCutTol = 1e-4;
        if (fabs(lb2-cut_lb) >= localCutTol || fabs(cut_ub-ub2) >= localCutTol) {
	  if (ub2 < infty) {
	    printf(" Strengthening ub of local cut for constraint %d from %e to %e\n", gindex, cut_ub, ub2);
	  }
	  else {
	    printf(" Strengthening ub of local cut for constraint %d from %e to %e\n", gindex, cut_lb, lb2);
	  }
	  // Now we generate a new cut
	  OsiRowCut newCut2;
	  newCut2.setEffectiveness(99.99e99);
	  newCut2.setLb(lb2);
	  newCut2.setUb(ub2);
	  newCut2.setRow(cut2);
	  cs.insert(newCut2);
        }
      }
    }
    return true;
  }

  bool CutStrengthener::StrengthenCut(SmartPtr<TMINLP> tminlp,
				      int constr_index,
				      const CoinPackedVector& row,
				      int n,
				      const double* x,
				      const double* x_l,
				      const double* x_u,
				      double& lb,
				      double& ub)
  {
    // Get data to set up the NLP to be solved
    Index nele_grad_gi;
    Index* jCol = new Index[n+1];
    bool new_x = true;
    if (constr_index == -1) {
      // Objective function
      // Compute random perturbation of point
      double* x_rand = new double[n];
      for (int i=0; i<n; i++) {
	const double radius = CoinMin(1., x_u[i]-x_l[i]);
	const double p = CoinMax(CoinMin(x[i]-0.5*radius, x_u[i]-radius),
				 x_l[i]);
	x_rand[i] = p + radius*CoinDrand48();
      }
      Number* grad_f = new Number[n];
      bool retval = tminlp->eval_grad_f(n, x_rand, new_x, grad_f);
      delete [] x_rand;
      if (!retval) {
	delete [] grad_f;
	delete [] jCol;
	return false;
      }
      nele_grad_gi = 0;
      for (int i=0; i<n; i++) {
	if (grad_f[i] != 0.) {
	  jCol[nele_grad_gi++] = i;
	}
      }
      delete [] grad_f;
      jCol[nele_grad_gi++] = n; // for the z variable
    }
    else {
      if (!tminlp->eval_grad_gi(n, x, new_x, constr_index, nele_grad_gi,
				jCol, NULL)) {
	delete [] jCol;
	return false;
      }
    }

    bool lower_bound;
    if (lb <= -COIN_DBL_MAX) {
      assert(ub < COIN_DBL_MAX);
      lower_bound = false;
    }
    else {
      assert(ub >= COIN_DBL_MAX);
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
    have_final_bound_(false),
    grad_f_(NULL)
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
	printf("There is an index (%d) in the cut that does not appear in the constraint.\n",idx);
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
      if (j < n) {
	x_l_[i] = x_l_orig[j];
	x_u_[i] = x_u_orig[j];
      }
      else {
	x_l_[i] = -1e100;
	x_u_[i] = 1e100;
      }
    }

    if (constr_index_ == -1) {
      grad_f_ = new Number[n_orig_];
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
    delete [] grad_f_;
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
    if (constr_index_ == -1) {
      g_l[0] = -1e100;
      g_u[0] = 0.;
    }
    else {
      Number* x_l_orig = new Number[n_orig_];
      Number* x_u_orig = new Number[n_orig_];
      Number* g_l_orig = new Number[m_orig_];
      Number* g_u_orig = new Number[m_orig_];

      if (!tminlp_->get_bounds_info(n_orig_, x_l_orig, x_u_orig,
				  m_orig_, g_l_orig, g_u_orig)) {
	delete [] x_l_orig;
	delete [] x_u_orig;
	delete [] g_l_orig;
	delete [] g_u_orig;
	return false;
      }

      g_l[0] = g_l_orig[constr_index_];
      g_u[0] = g_u_orig[constr_index_];

      delete [] x_l_orig;
      delete [] x_u_orig;
      delete [] g_l_orig;
      delete [] g_u_orig;
    }

    for (Index i=0; i<nvar_constr_; i++) {
      x_l[i] = x_l_[i];
      x_u[i] = x_u_[i];
    }

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
      if (constr_index_ == -1) {
	for (Index i=0; i<n-1; i++) {
	  x[i] = starting_point_[var_indices_[i]];
	}
	x[n-1] = 0.;
      }
      else {
	for (Index i=0; i<n; i++) {
	  x[i] = starting_point_[var_indices_[i]];
	}
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
    bool retval;
    if (constr_index_ == -1) {
      retval = tminlp_->eval_f(n_orig_, x_full_, new_x, g[0]);
      g[0] -= x[n-1];
    }
    else {
      retval = tminlp_->eval_gi(n_orig_, x_full_, new_x, constr_index_, g[0]);
    }
    return retval;
  }

  bool CutStrengthener::StrengtheningTNLP::
  eval_jac_g(Index n, const Number* x, bool new_x,
	     Index m, Index nele_jac, Index* iRow, Index *jCol,
	     Number* values)
  {
    bool retval = true;
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
      if (constr_index_ == -1) {
	retval = tminlp_->eval_grad_f(n_orig_, x_full_, new_x, grad_f_);
	if (retval) {
	  for (Index i=0; i<n-1; i++) {
	    values[i] = grad_f_[var_indices_[i]];
	  }
	  values[n-1] = -1.;
	}
      }
      else {
	retval = tminlp_->eval_grad_gi(n_orig_, x_full_, new_x, constr_index_,
				       nele_jac, NULL, values);
      }
    }
    return retval;
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
    if (constr_index_ == -1) {
      for (Index i=0; i<nvar_constr_-1; i++) {
	x_full_[var_indices_[i]] = x[i];
      }
    }
    else {
      for (Index i=0; i<nvar_constr_; i++) {
	x_full_[var_indices_[i]] = x[i];
      }
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
