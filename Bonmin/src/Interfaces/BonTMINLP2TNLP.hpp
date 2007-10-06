// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004

#ifndef __TMINLP2TNLP_HPP__
#define __TMINLP2TNLP_HPP__

#include "IpTNLP.hpp"
#include "BonTMINLP.hpp"
#include "IpSmartPtr.hpp"
#include "IpIpoptApplication.hpp"
#include "IpOptionsList.hpp"


namespace Bonmin
{
  class IpoptInteriorWarmStarter;
  /** This is an adapter class that converts a TMINLP to
   *  a TNLP to be solved by Ipopt. It allows an external
   *  caller to modify the bounds of variables, allowing
   *  the treatment of binary and integer variables as
   *  relaxed, or fixed
   */
  class TMINLP2TNLP : public Ipopt::TNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    TMINLP2TNLP(const SmartPtr<TMINLP> tminlp
#ifdef WARM_STARTER
        ,
        const OptionsList& options
#endif
        );

    /** Copy Constructor */
    TMINLP2TNLP(const TMINLP2TNLP&);

    /** Default destructor */
    virtual ~TMINLP2TNLP();
    //@}

    /** Copies the modification made to TNLP by the user (modifications
    such as changing bound changing starting point,...).
    this and other should be two instances of the same problem
    I am trying to mimic a copy construction for Cbc
    use with great care not safe.
    */
    void copyUserModification(const TMINLP2TNLP& other);


    /**@name Methods to modify the MINLP and form the NLP */
    //@{

    /** Get the number of variables */
    Index num_variables()
    {
      return n_;
    }

    /** Get the number of constraints */
    Index num_constraints()
    {
      return m_;
    }
    /** Get the nomber of nz in hessian */
    Index nnz_h_lag()
    {
      return nnz_h_lag_;
    }
    /** Get the variable types */
    const TMINLP::VariableType* var_types()
    {
      return var_types_;
    }

    //Print variable types to screen
    void printVarTypes()
    {
      //       for (int i = 0 ; i < n_ ; i++)
      // 	{
      // 	  std::cout<<i<<"\t"<<var_types_[i]<<std::endl;
      // 	}
    }

    /** Get the current values for the lower bounds */
    const Number* x_l()
    {
      return x_l_;
    }
    /** Get the current values for the upper bounds */
    const Number* x_u()
    {
      return x_u_;
    }

    /** Get the original values for the lower bounds */
    const Number* orig_x_l() const
    {
      return orig_x_l_;
    }
    /** Get the original values for the upper bounds */
    const Number* orig_x_u() const
    {
      return orig_x_u_;
    }

    /** Get the current values for constraints lower bounds */
    const Number* g_l()
    {
      return g_l_;
    }
    /** Get the current values for constraints upper bounds */
    const Number* g_u()
    {
      return g_u_;
    }

    /** get the starting primal point */
    const Number * x_init() const
    {
      return x_init_;
    }

    /** get the user provided starting primal point */
    const Number * x_init_user() const
    {
      return x_init_user_;
    }

    /** get the starting dual point */
    const Number * duals_init() const
    {
      return duals_init_;
    }

    /** get the solution values */
    const Number* x_sol() const
    {
      return x_sol_;
    }

    /** get the g solution (activities) */
    const Number* g_sol() const
    {
      return g_sol_;
    }

    /** get the dual values */
    const Number* duals_sol() const
    {
      return duals_sol_;
    }

    /** Get Optimization status */
    SolverReturn optimization_status() const
    {
      return return_status_;
    }

    /** Get the objective value */
    Number obj_value() const
    {
      return obj_value_;
    }

    /** Manually set objective value. */
    void set_obj_value(Number value)
    {
      obj_value_ = value;
    }

    /** force solution to be fractionnal.*/
    void force_fractionnal_sol();

    /** Change the bounds on the variables */
    void SetVariablesBounds(Index n,
                            const Number * x_l,
                            const Number * x_u);

    /** Change the lower bound on the variables */
    void SetVariablesLowerBounds(Index n,
                               const Number * x_l);

    /** Change the upper bound on the variable */
    void SetVariablesUpperBounds(Index n,
                                const Number * x_u);

    /** Change the bounds on the variable */
    void SetVariableBounds(Index var_no, Number x_l, Number x_u);

    /** Change the lower bound on the variable */
    void SetVariableLowerBound(Index var_no, Number x_l);

    /** Change the upper bound on the variable */
    void SetVariableUpperBound(Index var_no, Number x_u);

    /** Change the starting point */
    void SetStartingPoint(Index n, const Number* x_init);

    /** reset the starting point to original one. */
    void resetStartingPoint();

    /** Set component ind of the starting point. */
    void setxInit(Index ind,const Number val);

    /** set the starting point to x_init */
    void setxInit(Index n,const Number* x_init);

    /** Set component ind of the dual starting point. */
    void setDualInit(Index ind, const Number val);

    /** set the dual starting point to duals_init */
    void setDualsInit(Index n, const Number* duals_init);


    /** Set the contiuous solution */
    void Set_x_sol(Index n, const Number* x_sol);

    /** Change the type of the variable */
    void SetVariableType(Index n, TMINLP::VariableType type);
    //@}
    /** Procedure to ouptut relevant informations to reproduce a sub-problem.
      Compare the current problem to the problem to solve
      and writes files with bounds which have changed and current starting point.
      */
    void outputDiffs(const std::string& probName, const std::string* varNames);

    /**@name methods to gather information about the NLP */
    //@{
    /** This call is just passed onto the TMINLP object */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag,
        TNLP::IndexStyleEnum& index_style);

    /** The caller is allowed to modify the bounds, so this
     *  method returns the internal bounds information
     */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u);

    /** Returns the constraint linearity.
     * array should be alocated with length at least n. (default implementation
     *  just return false and does not fill the array).*/
    virtual bool get_constraints_linearity(Index m, LinearityType* const_types)
    {
      return tminlp_->get_constraints_linearity(m, const_types);
    }
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

    /** Method that returns scaling parameters. 
     */
    virtual bool get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling);


    /** Methat that returns an Ipopt IteratesVector that has the
     *  starting point for all internal varibles. */
    virtual bool get_warm_start_iterate(IteratesVector& warm_start_iterate);

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

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,
        Index n, const Number* x, const Number* z_L, const Number* z_U,
        Index m, const Number* g, const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq);
    /** Intermediate Callback method for the user.  Providing dummy
     *  default implementation.  For details see IntermediateCallBack
     *  in IpNLP.hpp. */
    virtual bool intermediate_callback(AlgorithmMode mode,
        Index iter, Number obj_value,
        Number inf_pr, Number inf_du,
        Number mu, Number d_norm,
        Number regularization_size,
        Number alpha_du, Number alpha_pr,
        Index ls_trials,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq);
    //@}

    /** Method called to check wether a problem has still some variable not fixed. If there are no more
        unfixed vars, checks wether the solution given by the bounds is feasible.*/

    /** @name Methods for setting and getting the warm starter */
    //@{
    void SetWarmStarter(SmartPtr<IpoptInteriorWarmStarter> warm_starter);

      SmartPtr<IpoptInteriorWarmStarter> GetWarmStarter();

    //@}
      
      /** Say if has a specific function to compute upper bounds*/
      virtual bool hasUpperBoundingObjective(){
        return tminlp_->hasUpperBoundingObjective();}

      /** Evaluate the upper bounding function at given point and store the result.*/
    void evaluateUpperBoundingFunction(const double * x);
    
    /** @name Methods for evaluation the linear cuts. */
    //@{
    /** Evaluate g for the linear cuts.
        \param n number of variables in the problem
        \param x points at which to evaluate
        \param nCuts number of cuts
        \param g where to put the results (array of size at least nCuts).
      */
    void eval_g_add_linear_cuts(Index n, const Number * x,
                                Index nCuts, Number* g);
    /** Evaluate the gradients of the cuts.*/
    void eval_jac_g_add_linear_cuts(Index nele_jac, Index* iRow,
				    Index *jCol, Number* values);
    //@}

  protected:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    TMINLP2TNLP();

    /** Overloaded Equals Operator */
    void operator=(const TMINLP2TNLP&);
    //@}

    /** pointer to the tminlp that is being adapted */
    SmartPtr<TMINLP> tminlp_;

    /** @name Internal copies of data allowing caller to modify the MINLP */
    //@{
    /// Number of variables
    Index n_;
    /// Number of constraints
    Index m_;
    /// Number of non-zeroes in the constraints jacobian.
    Index nnz_jac_g_;
    /// Number of non-zeroes in the lagrangian hessian
    Index nnz_h_lag_;
    /**index style (fortran or C)*/
    TNLP::IndexStyleEnum index_style_;
    /// Types of the variable (TMINLP::CONTINUOUS, TMINLP::INTEGER, TMINLP::BINARY).
    TMINLP::VariableType* var_types_;
    /// Current lower bounds on variables
    Number* x_l_;
    /// Current upper bounds on variables
    Number* x_u_;
    /// Original lower bounds on variables
    Number* orig_x_l_;
    /// Original upper bounds on variables
    Number* orig_x_u_;
    /// Lower bounds on constraints values
    Number* g_l_; // These are not modified, but a copy is kept for consistency
    /// Upper bounds on constraints values
    Number* g_u_;
    /// Initial primal point
    Number* x_init_;
    /** Initial values for all dual multipliers (constraints then lower bounds then upper bounds) */
    Number * duals_init_;
    /** storage capacity for starting point.*/
    Index capacity_x_init_;
    /// User-provideed initial prmal point
    Number* x_init_user_;
    /// Optimal solution
    Number* x_sol_;
    /// Activities of constraint g( x_sol_)
    Number * g_sol_;
    /** Dual multipliers of constraints and bounds*/
    Number* duals_sol_;

    /** Return status of the optimization process*/
    SolverReturn return_status_;
    /** Value of the optimal solution found by Ipopt */
    Number obj_value_;
    //@}

    /** @name Warmstart object and related data */
    //@{
    /** Pointer to object that holds warmstart information */
    SmartPtr<IpoptInteriorWarmStarter> curr_warm_starter_;
    /** Value for a lower bound that denotes -infinity */
    Number nlp_lower_bound_inf_;
    /** Value for a upper bound that denotes infinity */
    Number nlp_upper_bound_inf_;
    /** Option from Ipopt - we currently use it to see if we want to
     *  use some clever warm start or just the last iterate from the
     *  previous run */
    bool warm_start_entire_iterate_;
    /** Do we need a new warm starter object */
    bool need_new_warm_starter_;
    //@}


    /** Private method that throws an exception if the variable bounds
     * are not consistent with the variable type */
    void throw_exception_on_bad_variable_bound(Index i);

    private:
    void resizeStartingPoint();
  };

} // namespace Ipopt

#endif
