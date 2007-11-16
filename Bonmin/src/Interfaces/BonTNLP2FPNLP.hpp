// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
//
// Authors:  Pierre Bonami 06/10/2005

#ifndef _TNLP2FPNLP_HPP_
#define _TNLP2FPNLP_HPP_

#include "IpTNLP.hpp"
#include "BonTMINLP.hpp"
#include "IpSmartPtr.hpp"
#include "BonTypes.hpp"

namespace Bonmin
{
  /** This is an adapter class to convert an NLP to a Feasibility Pump NLP
   *  by changing the objective function to the (2-norm) distance to a point.
   * The extra function is set_dist2point_obj(int n, const double *, const int *)
   */
  class TNLP2FPNLP : public Ipopt::TNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    TNLP2FPNLP(const SmartPtr<TNLP> tnlp, double objectiveScalingFactor = 100);

    /** Default destructor */
    virtual ~TNLP2FPNLP();
    //@}

    /**@name Methods to change the objective function*/
    //@{
    /** \brief Set the point to which distance is minimized.
    * The distance is minimize in a subspace define by a subset of coordinates
     * \param n number of coordinates on which distance is minimized
     * \param inds indices of the coordinates on which distance is minimized
     * \param vals values of the point for coordinates in ind
     */
    void set_dist2point_obj(int n, const Number * vals, const Index * inds);

     /** Set the value for simgma */
     void setSigma(double sigma){
       assert(sigma >= 0.);
       sigma_ = sigma;}
     /** Set the value for lambda*/
     void setLambda(double lambda){
       assert(lambda >= 0. && lambda <= 1.);
       lambda_ = lambda;}
     /** Set the value for simgma */
     void setNorm(int norm){
       assert(norm >0 && norm < 3);
       norm_ = norm;}
    //@}

    /**@name methods to gather information about the NLP */
    //@{
    /** get info from nlp_ and add hessian information */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style);

    /** This call is just passed onto tnlp_
     */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u)
    {
      return tnlp_->get_bounds_info(n, x_l , x_u, m, g_l, g_u);
    }

    /** Passed onto tnlp_
     */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number* lambda)
    {
      return tnlp_->get_starting_point(n, init_x, x,
          init_z, z_L, z_U, m, init_lambda, lambda);
    }

    /** overloaded to return the value of the objective function */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
        Number& obj_value);

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
        Number* grad_f);

    /** passed onto tnlp_ */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
        Index m, Number* g)
    {
      return tnlp_->eval_g(n, x, new_x, m, g);
    }

    /** Passed onto tnlp_ */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values)
    {
      return tnlp_->eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }

    /** Evaluate the modified Hessian of the Lagrangian*/
    virtual bool eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,
        Index n, const Number* x, const Number* z_L, const Number* z_U,
        Index m, const Number* g, const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq);
    //@}

    /** @name Scaling of the objective function */
    //@{
    void setObjectiveScaling(double value)
    {
      objectiveScalingFactor_ = value;
    }
    double getObjectiveScaling() const
    {
      return objectiveScalingFactor_;
    }

  private:
    /** @name Internal methods to help compute the distance, its gradient and hessian */
    //@{
    /** Compute the norm-2 distance to the current point to which distance is minimized. */
    double dist2point(const Number *x);
    //@}
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    TNLP2FPNLP();

    /** Copy Constructor */
    TNLP2FPNLP(const TNLP2FPNLP&);

    /** Overloaded Equals Operator */
    void operator=(const TNLP2FPNLP&);
    //@}

    /** pointer to the tminlp that is being adapted */
    SmartPtr<TNLP> tnlp_;

    /** @name Data for storing the point the distance to which is minimized */
    //@{
    /// Indices of the variables for which distance is minimized (i.e. indices of integer variables in a feasibility pump setting)
    vector<Index> inds_;
    /// Values of the point to which we separate (if x is the point vals_[i] should be x[inds_[i]] )
    vector<Number> vals_;
    /** value for the convex combination to take between original objective and distance function.
      * ( take lambda_ * distance + (1-lambda) sigma f(x).*/
    double lambda_;
    /** Scaling for the original objective.*/
    double sigma_;
   /** Norm to use (L_1 or L_2).*/
   int norm_;
    //@}

    /// Scaling factor for the objective
    double objectiveScalingFactor_;

  };

} // namespace Ipopt

#endif /*_TNLP2FPNLP_HPP_*/
