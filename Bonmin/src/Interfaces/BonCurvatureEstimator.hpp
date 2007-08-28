// Copyright (C) 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:   Andreas Waechter                 IBM    2006-10-11

#ifndef __BONCURVATUREESTIMATOR_HPP__
#define __BONCURVATUREESTIMATOR_HPP__

#include "IpTNLP.hpp"
#include "IpOptionsList.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpCompoundVector.hpp"
#include <vector>

namespace Ipopt {
  //forward declarations
  class TSymLinearSolver;
}

namespace Bonmin
{
  using namespace Ipopt;

  /** 
   */
  class CurvatureEstimator: public ReferencedObject
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor.  It is given the options list to extract options
     *  specifying linear solver options. */
    CurvatureEstimator(
      SmartPtr<Journalist> jnlst,
      SmartPtr<OptionsList> options,
      SmartPtr<TNLP> tnlp);

    /** Destructor */
    virtual ~CurvatureEstimator();
    //@}

    /** Method for computing a direction projected_d related to the
	given direction orig_d and the two-sided product of
	projected_d with Hessian of Lagrangian.  The arrays x, y_c,
	and y_d constain the primal and dual variables definiting the
	Lagrangian Hessian. The vectors active_d and active_x contain
	the indices of active inequality and bound constraints,
	respectively.  A positive index number is interpreted to
	belong to an upper bound, and a negative number to a lower
	bound.  The return status is false if the computation was not
	possible, and true otherwise. */

    bool ComputeNullSpaceCurvature(
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
      Number& dTHLagd);

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    CurvatureEstimator();

    /** Copy Constructor */
    CurvatureEstimator(const CurvatureEstimator&);

    /** Overloaded Equals Operator */
    void operator=(const CurvatureEstimator&);
    //@}

    /** @name Items related to handling the linear solver object */
    //@{
    SmartPtr<Journalist> jnlst_;
    SmartPtr<OptionsList> options_;
    /** prefix to be used when parsion the options for the linear
     *  solver */
    std::string prefix_;
    /** Strategy object for solving the projection matrix for
	equality constraints only */
    SmartPtr<TSymLinearSolver> eq_tsymlinearsolver_;
    /** Strategy object for solving the projection matrix for
	all active constraints */
    SmartPtr<TSymLinearSolver> all_tsymlinearsolver_;
    //@}

    /** @name Information about the tnlp */
    //@{
    SmartPtr<TNLP> tnlp_;
    Index n_;
    Number* grad_f_;
    Index m_;
    Index nnz_jac_;
    Index* irows_jac_;
    Index* jcols_jac_;
    Number* jac_vals_;
    Index nnz_hess_;
    Index* irows_hess_;
    Index* jcols_hess_;
    Number* hess_vals_;
    //@}

    /** @name Information about activities for equality projection
	only */
    //@{
    /** Number of free x variables */
    Index eq_nx_free_;
    /** Map for pointing from the original x space to the one without
     *  fixed variables.  */
    Index* eq_x_free_map_;
    /** Number of active constraints */
    Index eq_ng_fixed_;
    /** Map for pointing from the original constraint space to the one
     *  with only active constraints */
    Index* eq_g_fixed_map_;
    //@}

    /** @name Information about activities for projection
	with respect to all active constraints */
    //@{
    /** Number of free x variables */
    Index all_nx_free_;
    /** Map for pointing from the original x space to the one without
     *  fixed variables.  */
    Index* all_x_free_map_;
    /** Number of active constraints */
    Index all_ng_fixed_;
    /** Map for pointing from the original constraint space to the one
     *  with only active constraints */
    Index* all_g_fixed_map_;
    //@}

    /** Space for most recent computed least-square multipliers */
    Number* lambda_;

    /** Space for storing the direction projected into the equality
     *  constraints only.  Having this array around all the time means
     *  that we don't have to be so careful about memory leaks. */
    Number* eq_projected_d_;

    /** @name Items for handling the projection system for equality
     *  constraints only */
    //@{
    /** Compound Matrix space for storing the linear system for the
     *  projection */
    SmartPtr<CompoundSymMatrixSpace> eq_comp_proj_matrix_space_;
    /** Compound Matrix storing the current projection matrix */
    SmartPtr<CompoundSymMatrix> eq_comp_proj_matrix_;
    /** Compound Vector space for storing right hand side and solution
     *  for the projection system */
    SmartPtr<CompoundVectorSpace> eq_comp_vec_space_;
    //@}

    /** @name Items for handling the projection system for all
     *  constraints */
    //@{
    /** Compound Matrix space for storing the linear system for the
     *  projection */
    SmartPtr<CompoundSymMatrixSpace> all_comp_proj_matrix_space_;
    /** Compound Matrix storing the current projection matrix */
    SmartPtr<CompoundSymMatrix> all_comp_proj_matrix_;
    /** Compound Vector space for storing right hand side and solution
     *  for the projection system */
    SmartPtr<CompoundVectorSpace> all_comp_vec_space_;
    //@}

    /** Storing the activities */
    //@{
    std::vector<int> active_x_;
    std::vector<int> active_g_;
    //@}

    bool initialized_;

    bool Initialize();

    bool PrepareNewMatrixStructure(
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
      SmartPtr<CompoundVectorSpace>& comp_vec_space);

    bool PrepareNewMatrixValues(
      const Index* x_free_map,
      const Index* g_fixed_map,
      SmartPtr<CompoundSymMatrixSpace>& comp_proj_matrix_space,
      SmartPtr<CompoundSymMatrix>& comp_proj_matrix,
      SmartPtr<TSymLinearSolver>& tsymlinearsolver);

    bool SolveSystem(
      const Number* rhs_x,
      const Number* rhs_g,
      Number* sol_x, Number* sol_g,
      const Index* x_free_map,
      const Index* g_fixed_map,
      SmartPtr<CompoundVectorSpace>& comp_vec_space,
      SmartPtr<CompoundSymMatrix>& comp_proj_matrix,
      SmartPtr<TSymLinearSolver>& tsymlinearsolver);

    bool Compute_dTHLagd(
      const Number* d, const Number* x, bool new_x, const Number* lambda,
      bool new_lambda,  Number& dTHLagd);
  };

} // namespace Ipopt
#endif
