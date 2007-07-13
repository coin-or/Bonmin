// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                   (derived from BonTMINLP2TNLP.hpp)            12/22/2006

#ifndef __BONBRANCHINGTQP_HPP__
#define __BONBRANCHINGTQP_HPP__

#include "BonTMINLP2TNLP.hpp"

namespace Bonmin
{
  /** This is an adapter class that converts a TMINLP2TNLP object into
   *  another such class, which is now just a QP.  The QP is the
   *  linear quadratic of the TNLP at the optimal point.  The purpose
   *  of the BranchingTQP is that it is used in a strong-branching
   *  framework, we strong branching is only done for the QP
   *  approximation of the TNLP, not on the TNLP itself.
   */
  class BranchingTQP : public TMINLP2TNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    BranchingTQP(const TMINLP2TNLP& tminlp2tnlp);

    /** Default destructor */
    virtual ~BranchingTQP();
    //@}

    /**@name methods to gather information about the NLP, only those
     *  that need to be overloaded from TMINLP2TNLP */
    //@{
    /** Method called by Ipopt to get the starting point. The bools
     *  init_x and init_lambda are both inputs and outputs. As inputs,
     *  they indicate whether or not the algorithm wants you to
     *  initialize x and lambda respectively. If, for some reason, the
     *  algorithm wants you to initialize these and you cannot, set
     *  the respective bool to false.
     */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number* lambda);

    /** Returns the value of the objective function in x*/
    virtual bool eval_f(Index n, const Number* x, bool new_x,
        Number& obj_value);

    /** Returns the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
        Number* grad_f);

    /** Returns the vector of constraint values in x*/
    virtual bool eval_g(Index n, const Number* x, bool new_x,
        Index m, Number* g);

    /** Returns the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values);

    /** Return the hessian of the
     *  lagrangian. The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  only */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values);
    //@}

    /** Accessor Methods for QP data */
    //@{
    const Number ObjVal()
    {
      return obj_val_;
    }
    const Number* ObjGrad()
    {
      return obj_grad_;
    }
    const Number* ObjHessVals()
    {
      return obj_hess_;
    }
    const Index* ObjHessIRow()
    {
      return obj_hess_irow_;
    }
    const Index* ObjHessJCol()
    {
      return obj_hess_jcol_;
    }
    const Number* ConstrRhs()
    {
      return g_vals_;
    }
    const Number* ConstrJacVals()
    {
      return g_jac_;
    }
    const Index* ConstrJacIRow()
    {
      return g_jac_irow_;
    }
    const Index* ConstrJacJCol()
    {
      return g_jac_jcol_;
    }
    //@}

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
    BranchingTQP();

    /** Copy Constructor */
    BranchingTQP(const BranchingTQP&);

    /** Overloaded Equals Operator */
    void operator=(const BranchingTQP&);
    //@}

    /** @name static information about the QP's constraints and
     *  objective function */
    //@{
    Number obj_val_;
    Number* obj_grad_;
    Number* obj_hess_;
    Index* obj_hess_irow_;
    Index* obj_hess_jcol_;
    Number* g_vals_;
    Number* g_jac_;
    Index* g_jac_irow_;
    Index* g_jac_jcol_;
    //@}

    /** Displacement with respect to x_sol_copy_ */
    Number* d_;

    /** Copy of original x_sol_.  x_sol_ is changed after the first QP
     *  has been solved once. */
    Number* x_sol_copy_;

    /** Copy of original duals_sol_.  duals_sol_ is changed after the
     *  first QP has been solved once. */
    Number* duals_sol_copy_;

    /** Method for updating the displacement */
    void update_displacement(const Number* x);
  };

} // namespace Ipopt

#endif
