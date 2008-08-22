// Copyright (C) 2006, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:   Andreas Waechter                 IBM    2006-10-11

#include "BonCurvatureEstimator.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpGenTMatrix.hpp"
#include "IpIdentityMatrix.hpp"
#include "IpZeroMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"
#include <cmath>
#ifdef HAVE_MA27
# include "IpMa27TSolverInterface.hpp"
#endif
#ifdef HAVE_MA57
# include "IpMa57TSolverInterface.hpp"
#endif
#ifdef HAVE_MC19
# include "IpMc19TSymScalingMethod.hpp"
#endif
#ifdef HAVE_PARDISO
# include "IpPardisoSolverInterface.hpp"
#endif
#ifdef HAVE_TAUCS
# include "IpTAUCSSolverInterface.hpp"
#endif
#ifdef HAVE_WSMP
# include "IpWsmpSolverInterface.hpp"
#endif
#ifdef HAVE_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif

namespace Bonmin
{
  using namespace Ipopt;

  // ToDo: Consider NLP scaling?

  CurvatureEstimator::CurvatureEstimator
    (SmartPtr<Journalist> jnlst,
     SmartPtr<OptionsList> options,
     SmartPtr<TNLP> tnlp)
      :
      jnlst_(jnlst),
      options_(options),
      prefix_(""),
      tnlp_(tnlp),
      grad_f_(NULL),
      irows_jac_(NULL),
      jcols_jac_(NULL),
      jac_vals_(NULL),
      irows_hess_(NULL),
      jcols_hess_(NULL),
      hess_vals_(NULL),
      eq_x_free_map_(NULL),
      eq_g_fixed_map_(NULL),
      all_x_free_map_(NULL),
      all_g_fixed_map_(NULL),
      lambda_(NULL),
      eq_projected_d_(NULL),
      initialized_(false)
  {
    DBG_ASSERT(IsValid(jnlst));
    DBG_ASSERT(IsValid(options));
    DBG_ASSERT(IsValid(tnlp));

    ////////////////////////////////////////////////////////////////////
    // Create a strategy object for solving the linear system for the //
    // projection matrix                                              //
    ////////////////////////////////////////////////////////////////////

    // The following linear are taken from AlgBuilder in Ipopt
    SmartPtr<SparseSymLinearSolverInterface> SolverInterface1;
    SmartPtr<SparseSymLinearSolverInterface> SolverInterface2;
    std::string linear_solver;
    options->GetStringValue("linear_solver", linear_solver, prefix_);
    if (linear_solver=="ma27") {
#ifdef HAVE_MA27
      SolverInterface1 = new Ma27TSolverInterface();
      SolverInterface2 = new Ma27TSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MA27 not available.");
#endif

    }
    else if (linear_solver=="ma57") {
#ifdef HAVE_MA57
      SolverInterface1 = new Ma57TSolverInterface();
      SolverInterface2 = new Ma57TSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MA57 not available.");
#endif

    }
    else if (linear_solver=="pardiso") {
#ifdef HAVE_PARDISO
      SolverInterface1 = new PardisoSolverInterface();
      SolverInterface2 = new PardisoSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver Pardiso not available.");
#endif

    }
    else if (linear_solver=="taucs") {
#ifdef HAVE_TAUCS
      SolverInterface1 = new TAUCSSolverInterface();
      SolverInterface2 = new TAUCSSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver TAUCS not available.");
#endif

    }
    else if (linear_solver=="wsmp") {
#ifdef HAVE_WSMP
      SolverInterface1 = new WsmpSolverInterface();
      SolverInterface2 = new WsmpSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver WSMP not available.");
#endif

    }
    else if (linear_solver=="mumps") {
#ifdef HAVE_MUMPS
      SolverInterface1 = new MumpsSolverInterface();
      SolverInterface2 = new MumpsSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MUMPS not available.");
#endif

    }

    SmartPtr<TSymScalingMethod> ScalingMethod1;
    SmartPtr<TSymScalingMethod> ScalingMethod2;
    std::string linear_system_scaling;
    if (!options->GetStringValue("linear_system_scaling",
				 linear_system_scaling, prefix_)) {
      // By default, don't use mc19 for non-HSL solvers
      if (linear_solver!="ma27" && linear_solver!="ma57") {
        linear_system_scaling="none";
      }
    }
    if (linear_system_scaling=="mc19") {
#ifdef HAVE_MC19
      ScalingMethod1 = new Mc19TSymScalingMethod();
      ScalingMethod2 = new Mc19TSymScalingMethod();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear system scaling method MC19 not available.");
#endif

    }

    eq_tsymlinearsolver_ =
      new TSymLinearSolver(SolverInterface1, ScalingMethod1);
    all_tsymlinearsolver_ =
      new TSymLinearSolver(SolverInterface2, ScalingMethod2);
    // End of lines from AlgBuilder
  }

  CurvatureEstimator::~CurvatureEstimator()
  {
    if (initialized_) {
      delete [] grad_f_;
      delete [] irows_jac_;
      delete [] jcols_jac_;
      delete [] jac_vals_;
      delete [] irows_hess_;
      delete [] jcols_hess_;
      delete [] hess_vals_;
      delete [] eq_x_free_map_;
      delete [] eq_g_fixed_map_;
      delete [] all_x_free_map_;
      delete [] all_g_fixed_map_;
      delete [] eq_projected_d_;
      delete [] lambda_;
    }
  }

  bool CurvatureEstimator::Initialize()
  {
    DBG_ASSERT(!initialized_);

    //////////////////////////////////////
    // Prepare internal data structures //
    //////////////////////////////////////

    // Get sizes
    TNLP::IndexStyleEnum index_style;
    if (!tnlp_->get_nlp_info(n_, m_, nnz_jac_, nnz_hess_, index_style)) {
      return false;
    }

    // Space for gradient
    delete [] grad_f_;
    grad_f_ = NULL;
    grad_f_ = new Number[n_];

    // Get nonzero entries in the matrices
    delete [] irows_jac_;
    delete [] jcols_jac_;
    irows_jac_ = NULL;
    jcols_jac_ = NULL;
    irows_jac_ = new Index[nnz_jac_];
    jcols_jac_ = new Index[nnz_jac_];
    if (!tnlp_->eval_jac_g(n_, NULL, false, m_, nnz_jac_,
			   irows_jac_, jcols_jac_, NULL)) {
      return false;
    }
    if (index_style == TNLP::FORTRAN_STYLE) {
      for (Index i=0; i<nnz_jac_; i++) {
	irows_jac_[i]--;
	jcols_jac_[i]--;
      }
    }
    delete [] jac_vals_;
    jac_vals_ = NULL;
    jac_vals_ = new Number[nnz_jac_];

    delete [] irows_hess_;
    irows_hess_ = NULL;
    irows_hess_ = new Index[nnz_hess_];
    delete [] jcols_hess_;
    jcols_hess_ = NULL;
    jcols_hess_ = new Index[nnz_hess_];
    if (!tnlp_->eval_h(n_, NULL, false, 1., m_, NULL, false, nnz_hess_,
		       irows_hess_, jcols_hess_, NULL)) {
      return false;
    }
    if (index_style == TNLP::FORTRAN_STYLE) {
      for (Index i=0; i<nnz_hess_; i++) {
	irows_hess_[i]--;
	jcols_hess_[i]--;
      }
    }
    delete [] hess_vals_;
    hess_vals_ = NULL; // We set it to NULL, so that we know later
    // that we still need to compute the values

    // Get space for the activities maps
    delete [] eq_x_free_map_;
    delete [] eq_g_fixed_map_;
    delete [] all_x_free_map_;
    delete [] all_g_fixed_map_;
    eq_x_free_map_ = NULL;
    eq_g_fixed_map_ = NULL;
    all_x_free_map_ = NULL;
    all_g_fixed_map_ = NULL;
    eq_x_free_map_ = new Index[n_];
    eq_g_fixed_map_ = new Index[m_];
    all_x_free_map_ = new Index[n_];
    all_g_fixed_map_ = new Index[m_];

    // Get space for the multipliers
    delete [] lambda_;
    lambda_ = NULL;
    lambda_ = new Number[m_];

    // Get space for projected d
    delete [] eq_projected_d_;
    eq_projected_d_ = NULL;
    eq_projected_d_ = new Number[n_];

    initialized_ = true;
    return true;
  }

  bool
  CurvatureEstimator::ComputeNullSpaceCurvature(
    int n,
    const Number* x,
    bool new_x,
    const Number* x_l,
    const Number* x_u,
    const Number* g_l,
    const Number* g_u,
    bool new_bounds,
    const Number* z_L,
    const Number* z_U,
    int m,
    const Number* lam,
    bool new_mults,
    const Number* orig_d,
    Number* projected_d,
    Number& gradLagTd,
    Number& dTHLagd)
  {
#if 0
    printf("new_bounds = %d new_x = %d new_mults = %d\n", new_bounds, new_x, new_mults);
  for (int i=0;  i<n; i++) {
    printf("x[%3d] = %15.8e orig_d[%3d] = %15.8e z_L[%3d] = %15.8e z_U[%3d] = %15.8e\n",i,x[i],i,orig_d[i],i,z_L[i],i,z_U[i]);
  }
  for (int i=0; i<m; i++) {
    printf("lam[%3d] = %15.8e\n", i, lam[i]);
  }
#endif
    if (!initialized_) {
      Initialize();
      new_bounds = true;
      new_mults = true;
    }
    //DELETEME
    new_bounds = true;

    DBG_ASSERT(n == n_);
    DBG_ASSERT(m == m_);

    // If necessary, get new Jacobian values (for the original matrix)
    if (new_x) {
      if (!tnlp_->eval_jac_g(n_, x, new_x, m_, nnz_jac_,
			     NULL, NULL, jac_vals_)) {
	return false;
      }
    }

    // First we compute the direction projected into the space of only
    // the equality constraints
    if (new_x) {
      std::vector<int> dummy_active_x;
      std::vector<int> dummy_active_g;
      if (!PrepareNewMatrixStructure(x_l, x_u, g_l, g_u,
				     dummy_active_x, dummy_active_g,
				     eq_nx_free_, eq_x_free_map_,
				     eq_ng_fixed_, eq_g_fixed_map_,
				     eq_comp_proj_matrix_space_,
				     eq_comp_vec_space_)) {
	return false;
      }

      if (!PrepareNewMatrixValues(eq_x_free_map_, eq_g_fixed_map_,
				  eq_comp_proj_matrix_space_,
				  eq_comp_proj_matrix_,
				  eq_tsymlinearsolver_)) {
	return false;
      }
    }

    if (eq_ng_fixed_>0) {
      // Compute the projection of the direction
      if (!SolveSystem(orig_d, NULL, eq_projected_d_, NULL,
		       eq_x_free_map_, eq_g_fixed_map_,
		       eq_comp_vec_space_, eq_comp_proj_matrix_,
		       eq_tsymlinearsolver_)) {
	return false;
      }
      orig_d = eq_projected_d_; // This way we don't need to rememeber
				// how the direction was projected;
    }

    // If necessary, determine new activities
    bool new_activities = false;
    if (new_bounds || new_mults) {
      new_activities = true;
      active_x_.clear();
      if (eq_ng_fixed_>0) {
	const Number zTol = 1e-4;
	jnlst_->Printf(J_MOREDETAILED, J_NLP,
		       "List of variables considered fixed (with orig_d and z)\n");
	for (Index i=0; i<n; i++) {
	  if (x_l[i] < x_u[i]) {
	    if (orig_d[i]>0. && z_U[i]*orig_d[i]>zTol) {
	      active_x_.push_back(i+1);
	      jnlst_->Printf(J_MOREDETAILED, J_NLP,
			     "x[%5d] (%e,%e)\n", i, orig_d[i], z_U[i]);
	      DBG_ASSERT(x_u[i] < 1e19);
	    }
	    else if (orig_d[i]<0. && -z_L[i]*orig_d[i]>zTol) {
	      active_x_.push_back(-(i+1));
	      jnlst_->Printf(J_MOREDETAILED, J_NLP,
			     "x[%5d] (%e,%e)\n", i, orig_d[i], z_L[i]);
	      DBG_ASSERT(x_l[i] > -1e19);
	    }
	  }
	}
      }

      active_g_.clear();
      // Compute the product of the direction with the constraint Jacobian
      // This could be done more efficient if we have d in sparse format
      Number* jacTd = new Number[m];
      const Number zero = 0.;
      IpBlasDcopy(m, &zero, 0, jacTd, 1);
      for (Index i=0; i<nnz_jac_; i++) {
	const Index& irow = irows_jac_[i];
	const Index& jcol = jcols_jac_[i];
	jacTd[irow] += jac_vals_[i]*orig_d[jcol];
      }

      const Number lamTol = 1e-4;
      jnlst_->Printf(J_MOREDETAILED, J_NLP,
		     "List of constraints considered fixed (with lam and jacTd)\n");
      for (Index j=0; j<m; j++) {
	if (g_l[j] < g_u[j] && fabs(lam[j]) > lamTol) {
	  if (lam[j]*jacTd[j] > 0) {
	    if (lam[j] < 0.) {
	      active_g_.push_back(-(j+1));
	      DBG_ASSERT(g_l[j] > -1e19);
	    }
	    else {
	      active_g_.push_back(j+1);
	      DBG_ASSERT(g_u[j] < 1e19);
	    }
	    //	    active_g_.push_back(j+1);
	    jnlst_->Printf(J_MOREDETAILED, J_NLP,
			   "g[%5d] (%e,%e)\n", j, lam[j], jacTd[j]);
	  }
	}
      }
      delete [] jacTd;
    }

    // Check if the structure of the matrix has changed
    if (new_activities) {
      if (!PrepareNewMatrixStructure(x_l, x_u, g_l, g_u,
				     active_x_, active_g_,
				     all_nx_free_, all_x_free_map_,
				     all_ng_fixed_, all_g_fixed_map_,
				     all_comp_proj_matrix_space_,
				     all_comp_vec_space_)) {
	return false;
      }
    }

    bool new_lambda = false;
    if (new_x || new_activities) {
      if (!PrepareNewMatrixValues(all_x_free_map_, all_g_fixed_map_,
				  all_comp_proj_matrix_space_,
				  all_comp_proj_matrix_,
				  all_tsymlinearsolver_)) {
	return false;
      }

#ifdef lambdas
      // Compute least square multipliers for the given activities
      if (!tnlp_->eval_grad_f(n_, x, new_x, grad_f_)) {
	return false;
      }
      if (!SolveSystem(grad_f_, NULL, NULL, lambda_)) {
	return false;
      }
      IpBlasDscal(m_, -1., lambda_, 1);
      if (jnlst_->ProduceOutput(J_MOREVECTOR, J_NLP)) {
	jnlst_->Printf(J_MOREVECTOR, J_NLP,
		       "Curvature Estimator: least square multiplier:\n");
	for (Index i=0; i<m_; i++) {
	  jnlst_->Printf(J_MOREVECTOR, J_NLP, "lambda[%5d] = %23.16e\n",
			 i, lambda_[i]);
	}
      }
      new_lambda = true;
#endif
    }

    // Compute the projection of the direction
    if (!SolveSystem(orig_d, NULL, projected_d, NULL,
		     all_x_free_map_, all_g_fixed_map_,
		     all_comp_vec_space_, all_comp_proj_matrix_,
		     all_tsymlinearsolver_)) {
      return false;
    }

    // Sanity check to see if the thing is indeed in the null space
    // (if the constraint gradients are rank-deficient, the solver
    // might not have done a good job)
    Number* trash = new Number[m_];
    for (Index j=0; j<m_; j++) {
      trash[j] = 0.;
    }
    for (Index i=0; i<nnz_jac_; i++) {
      const Index &irow = irows_jac_[i];
      const Index &jcol = jcols_jac_[i];
      if (all_x_free_map_[jcol] >= 0 && all_g_fixed_map_[irow] >= 0) {
	trash[irow] += jac_vals_[i]*projected_d[jcol];
      }
    }
    if (jnlst_->ProduceOutput(J_MOREVECTOR, J_NLP)) {    
      for (Index j=0; j<m_; j++) {
	jnlst_->Printf(J_MOREVECTOR, J_NLP,
		       "nullspacevector[%5d] = %e\n", j, trash[j]);
      }
    }
    Index imax = IpBlasIdamax(m_, trash, 1)-1;
    Number max_trash = trash[imax];
    delete [] trash;
    const Number max_trash_tol = 1e-6;
    if (max_trash > max_trash_tol) {
      jnlst_->Printf(J_WARNING, J_NLP,
		     "Curvature Estimator: Bad solution from linear solver with max_red = %e:\n", max_trash);
      return false;
    }

    if (jnlst_->ProduceOutput(J_MOREVECTOR, J_NLP)) {
      jnlst_->Printf(J_MOREVECTOR, J_NLP,
		     "Curvature Estimator: original and projected directions are:\n");
      for (Index i=0; i<n_; i++) {
	jnlst_->Printf(J_MOREVECTOR, J_NLP,
		       "orig_d[%5d] = %23.16e proj_d[%5d] = %23.16e\n",
		       i, orig_d[i], i, projected_d[i]);
      }
    }

    gradLagTd = 0.;
#ifdef lambdas
    // Compute the product with the gradient of the Lagrangian
    gradLagTd = IpBlasDdot(n, projected_d, 1, grad_f_, 1);
    for (Index i=0; i<nnz_jac_; i++) {
      const Index &irow = irows_jac_[i];
      const Index &jcol = jcols_jac_[i];
      gradLagTd += lambda_[irow]*jac_vals_[i]*projected_d[jcol];
    }
#endif

    // Compute the product with the Hessian of the Lagrangian
    //    if (!Compute_dTHLagd(projected_d, x, new_x, lambda_, new_lambda, dTHLagd)) {
    if (!Compute_dTHLagd(projected_d, x, new_x, lam, new_lambda, dTHLagd)) {
      return false;
    }

#if 0
    printf("gradLagTd = %e dTHLagd = %e\n",gradLagTd,dTHLagd);
#endif
    return true;
  }

  bool CurvatureEstimator::PrepareNewMatrixStructure(
    const Number* x_l,
    const Number* x_u,
    const Number* g_l,
    const Number* g_u,
    std::vector<int>& active_x,
    std::vector<int>& active_g,
    Index& nx_free,
    Index* x_free_map,
    Index& ng_fixed,
    Index* g_fixed_map,
    SmartPtr<CompoundSymMatrixSpace>& comp_proj_matrix_space,
    SmartPtr<CompoundVectorSpace>& comp_vec_space)
  {
    // Deterimine which variables are free
    for (Index i=0; i<n_; i++) {
      x_free_map[i] = 0;
    }
    // fixed by activities
    for (std::vector<int>::iterator i=active_x.begin();
	 i != active_x.end(); i++) {
      DBG_ASSERT(*i != 0 && *i<=n_ && *i>=-n_);
      if (*i<0) {
	x_free_map[(-*i)-1] = -1;
	DBG_ASSERT(x_l[(-*i)-1] > -1e19);
      }
      else {
	x_free_map[(*i)-1] = -1;
	DBG_ASSERT(x_u[(*i)-1] < 1e19);
      }
    }
    // fixed by bounds
    for (Index i=0; i<n_; i++) {
      if (x_l[i] == x_u[i]) {
	x_free_map[i] = -1;
      }
    }
    // Correct the numbering in the x map and determine number of
    // free variables
    nx_free = 0;
    for (Index i=0; i<n_; i++) {
      if (x_free_map[i] == 0) {
	x_free_map[i] = nx_free++;
      }
    }

    // Determine which constraints are fixed
    for (Index j=0; j<m_; j++) {
      g_fixed_map[j] = -1;
    }
    // fixed by activities
    for (std::vector<int>::iterator i=active_g.begin();
	 i != active_g.end(); i++) {
      DBG_ASSERT(*i != 0 && *i<=m_ && *i>=-m_);
      if (*i<0) {
	g_fixed_map[(-*i)-1] = 0;
	DBG_ASSERT(g_l[(-*i)-1] > -1e19); //ToDo look at option?
      }
      else {
	g_fixed_map[(*i)-1] = 0;
	DBG_ASSERT(g_u[(*i)-1] < 1e19);
      }
    }
    // fixed by bounds
    for (Index j=0; j<m_; j++) {
      if (g_l[j] == g_u[j]) {
	g_fixed_map[j] = 0;
      }
    }
    // Correct the numbering in the g map and determine number of
    // fixed constraints
    ng_fixed = 0;
    for (Index j=0; j<m_; j++) {
      if (g_fixed_map[j] == 0) {
	g_fixed_map[j] = ng_fixed++;
      }
    }

    // Determine the row and column indices for the Jacobian of the fixed
    // constraints in the space of the free variables
    Index* iRows = new Index[nnz_jac_];
    Index* jCols = new Index[nnz_jac_];
    Index nnz_proj_jac = 0;
    for (Index i=0; i<nnz_jac_; i++) {
      const Index &irow = irows_jac_[i];
      const Index &jcol = jcols_jac_[i];
      if (x_free_map[jcol] >= 0 && g_fixed_map[irow] >= 0) {
	iRows[nnz_proj_jac] = g_fixed_map[irow] + 1;
	jCols[nnz_proj_jac] = x_free_map[jcol] + 1;
	nnz_proj_jac++;
      }
    }

    // Create the matrix space for the Jacobian matrices
    SmartPtr<GenTMatrixSpace> proj_jac_space =
      new GenTMatrixSpace(ng_fixed, nx_free, nnz_proj_jac, iRows, jCols);
    delete [] iRows;
    delete [] jCols;

    // Create Matrix space for the projection matrix
    comp_proj_matrix_space =
      new CompoundSymMatrixSpace(2, nx_free+ng_fixed);
    comp_proj_matrix_space->SetBlockDim(0, nx_free);
    comp_proj_matrix_space->SetBlockDim(1, ng_fixed);
    SmartPtr<SymMatrixSpace> identity_space =
      new IdentityMatrixSpace(nx_free);
    comp_proj_matrix_space->SetCompSpace(0, 0, *identity_space, true);
    comp_proj_matrix_space->SetCompSpace(1, 0, *proj_jac_space, true);

    // Create a Vector space for the rhs and sol
    comp_vec_space = new CompoundVectorSpace(2, nx_free+ng_fixed);
    SmartPtr<DenseVectorSpace> x_space = new DenseVectorSpace(nx_free);
    comp_vec_space->SetCompSpace(0, *x_space);
    SmartPtr<DenseVectorSpace> g_space = new DenseVectorSpace(ng_fixed);
    comp_vec_space->SetCompSpace(1, *g_space);

    return true;
  }

  bool CurvatureEstimator::PrepareNewMatrixValues(
    const Index* x_free_map,
    const Index* g_fixed_map,
    SmartPtr<CompoundSymMatrixSpace>& comp_proj_matrix_space,
    SmartPtr<CompoundSymMatrix>& comp_proj_matrix,
    SmartPtr<TSymLinearSolver>& tsymlinearsolver)
  {
    comp_proj_matrix = comp_proj_matrix_space->MakeNewCompoundSymMatrix();
    SmartPtr<Matrix> jac = comp_proj_matrix->GetCompNonConst(1, 0);
    SmartPtr<GenTMatrix> tjac = static_cast<GenTMatrix*> (GetRawPtr(jac));
    Number* vals = tjac->Values();
    Index inz=0;
    for (Index i=0; i<nnz_jac_; i++) {
      Index irow = irows_jac_[i];
      Index jcol = jcols_jac_[i];
      if (x_free_map[jcol] >= 0 && g_fixed_map[irow] >= 0) {
	vals[inz++] = jac_vals_[i];
      }
    }
    DBG_ASSERT(inz == tjac->Nonzeros());

    // We need to reset the linear solver object, so that it knows
    // that the structure of the linear system has changed
    tsymlinearsolver->ReducedInitialize(*jnlst_, *options_, prefix_);

    return true;
  }

  bool CurvatureEstimator::SolveSystem(
    const Number* rhs_x,
    const Number* rhs_g,
    Number* sol_x, Number* sol_g,
    const Index* x_free_map,
    const Index* g_fixed_map,
    SmartPtr<CompoundVectorSpace>& comp_vec_space,
    SmartPtr<CompoundSymMatrix>& comp_proj_matrix,
    SmartPtr<TSymLinearSolver>& tsymlinearsolver)
  {
    // Create a vector for the RHS
    SmartPtr<CompoundVector> rhs = comp_vec_space->MakeNewCompoundVector();
    SmartPtr<Vector> vrhs_x = rhs->GetCompNonConst(0);
    SmartPtr<Vector> vrhs_g = rhs->GetCompNonConst(1);
    // Now fill this vector with the values, extracting the relevant entries
    if (rhs_x) {
      SmartPtr<DenseVector> drhs_x =
	static_cast<DenseVector*> (GetRawPtr(vrhs_x));
      Number* xvals = drhs_x->Values();
      for (Index i=0; i<n_; i++) {
	const Index& ix = x_free_map[i];
	if (ix>=0) {
	  xvals[ix] = rhs_x[i];
	}
      }
    }
    else {
      vrhs_x->Set(0.);
    }
    if (rhs_g) {
      SmartPtr<DenseVector> drhs_g =
	static_cast<DenseVector*> (GetRawPtr(vrhs_g));
      Number* gvals = drhs_g->Values();
      for (Index j=0; j<m_; j++) {
	const Index& jg = g_fixed_map[j];
	if (jg>=0) {
	  gvals[jg] = rhs_g[j];
	}
      }
    }
    else {
      vrhs_g->Set(0.);
    }

    // Solve the linear system
    SmartPtr<CompoundVector> sol = comp_vec_space->MakeNewCompoundVector();
    ESymSolverStatus solver_retval =
      tsymlinearsolver->Solve(*comp_proj_matrix, *rhs, *sol, false, 0);
    // For now we try, if the system is reported to be singular, to
    // still get a solution; ToDo
    if (solver_retval == SYMSOLVER_SINGULAR) {
      jnlst_->Printf(J_DETAILED, J_NLP,
		     "Projection matrix reported to be singular, but we try to obtain a solution anyway.\n");
      solver_retval =
	tsymlinearsolver->Solve(*comp_proj_matrix, *rhs, *sol, false, 0);
    }

    if (solver_retval != SYMSOLVER_SUCCESS) {
      // DELETEME
      //printf("Return code from Solver is %d\n", solver_retval);
      return false;
    }

    // and get the solution out of it
    if (sol_x) {
      SmartPtr<Vector> vsol_x = sol->GetCompNonConst(0);
      SmartPtr<const DenseVector> dsol_x =
	static_cast<const DenseVector*> (GetRawPtr(vsol_x));
      const Number* xvals = dsol_x->Values();
      for (Index i=0; i<n_; i++) {
	const Index& ix = x_free_map[i];
	if (ix >=0) {
	  sol_x[i] = xvals[ix];
	}
	else {
	  sol_x[i] = 0.;
	}
      }
    }
    if (sol_g) {
      SmartPtr<Vector> vsol_g = sol->GetCompNonConst(1);
      SmartPtr<const DenseVector> dsol_g =
	static_cast<const DenseVector*> (GetRawPtr(vsol_g));
      const Number* gvals = dsol_g->Values();
      for (Index j=0; j<m_; j++) {
	const Index& ig = g_fixed_map[j];
	if (ig>=0) {
	  sol_g[j] = gvals[ig];
	}
	else {
	  sol_g[j] = 0.;
	}
      }
    }

    return true;
  }

  bool CurvatureEstimator::Compute_dTHLagd(
    const Number* d, const Number* x, bool new_x, const Number* lambda,
    bool new_lambda,  Number& dTHLagd)
  {
    if (new_x || new_lambda || !hess_vals_) {
      if (!hess_vals_) {
	hess_vals_ = new Number[nnz_hess_];
      }
      if (!tnlp_->eval_h(n_, x, new_x, 1., m_, lambda, new_lambda, nnz_hess_,
			 NULL, NULL, hess_vals_)) {
	return false;
      }
    }
    dTHLagd = 0.;
    for (Index i=0; i<nnz_hess_; i++) {
      Index irow = irows_hess_[i];
      Index jcol = jcols_hess_[i];
      if (irow == jcol) {
	dTHLagd += d[irow]*d[irow]*hess_vals_[i];
      }
      else {
	dTHLagd += 2.*d[irow]*d[jcol]*hess_vals_[i];
      }
    }
    return true;
  }


} // namespace Bonmin
